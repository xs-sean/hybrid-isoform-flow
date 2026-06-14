from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .config import ConfigError, Sample


IDP_OUTPUTS = ["isoform_detection.gpd", "isoform_prediction.gpd", "isoform.gpd"]
DEFAULT_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def make_idp_configs(
    samples: list[Sample],
    references: dict[str, Any],
    idp: dict[str, Any],
    out_dir: Path,
    sample_id: str | None = None,
    chromosome: str | None = None,
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for sample in samples:
        if sample_id and sample.sample_id != sample_id:
            continue
        suffix = f".{chromosome}" if chromosome else ""
        run_name = f"{sample.sample_id}{suffix}"
        values = {
            **idp,
            "temp_foldername": f"work/idp/{run_name}/temp",
            "output_foldername": f"results/idp/{run_name}",
            "genome_pathfilename": references["genome_fasta"],
            "gmap_index_pathfoldername": references.get("gmap_database", references.get("gmap_index", "")),
            "ref_annotation_pathfilename": references["reference_gpd"],
            "allref_annotation_pathfilename": references["allref_gpd"],
            "SR_jun_pathfilename": _sample_junction_bed(sample, chromosome),
            "SR_sam_pathfilename": _sample_sr_sam(sample, chromosome),
        }
        long_read = _sample_long_read(sample, chromosome)
        long_read_key = _long_read_config_key(long_read)
        values[long_read_key] = long_read
        config_path = out_dir / f"{run_name}.run.cfg"
        config_path.write_text(_format_idp_config(values))
        written.append(config_path)
    if sample_id and not written:
        raise ConfigError(f"No sample with sample_id {sample_id!r}")
    return written


def run_idp_configs(
    config_dir: Path,
    idp_root: Path,
    container_image: str,
    workdir: Path,
    dry_run: bool,
    config_file: Path | None = None,
) -> list[list[str]]:
    configs = [config_file] if config_file else sorted(config_dir.glob("*.run.cfg"))
    if not configs:
        raise ConfigError(f"No *.run.cfg files found in {config_dir}")
    run_idp = idp_root / "bin" / "runIDP.py"
    if not run_idp.exists():
        raise ConfigError(f"Missing upstream runIDP.py at {run_idp}")
    commands: list[list[str]] = []
    for cfg in configs:
        if container_image:
            run_idp_container = _container_path(run_idp, workdir)
            cfg_container = _container_path(cfg, workdir)
            command = [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{workdir.resolve()}:/work",
                "-w",
                "/work",
                container_image,
                "python2",
                run_idp_container,
                cfg_container,
                "0",
            ]
        else:
            command = ["python2", str(run_idp), str(cfg), "0"]
        commands.append(command)
        if not dry_run:
            subprocess.run(command, check=True)
    return commands


def _container_path(path: Path, workdir: Path) -> str:
    resolved_path = path.resolve()
    resolved_workdir = workdir.resolve()
    try:
        return str(Path("/work") / resolved_path.relative_to(resolved_workdir))
    except ValueError:
        return str(resolved_path)


def merge_gpd_outputs(idp_output_dir: Path, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for output_name in IDP_OUTPUTS:
        target = out_dir / output_name
        with target.open("w") as out:
            for source in sorted(idp_output_dir.glob(f"*/{output_name}")):
                text = source.read_text()
                out.write(text)
                if not text.endswith("\n"):
                    out.write("\n")
        written.append(target)
    return written


def merge_sample_gpd_outputs(idp_output_dir: Path, out_dir: Path, samples: list[Sample], chromosomes: list[str]) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for sample in samples:
        for output_name in IDP_OUTPUTS:
            target = out_dir / f"{sample.sample_id}.{output_name}"
            with target.open("w") as out:
                for chrom in chromosomes:
                    source = idp_output_dir / f"{sample.sample_id}.{chrom}" / output_name
                    if not source.exists():
                        continue
                    text = source.read_text()
                    out.write(text)
                    if not text.endswith("\n"):
                        out.write("\n")
            written.append(target)
    return written


def gpd_to_gtf(gpd_dir: Path, out_dir: Path, gene_pred_to_gtf: Path, dry_run: bool) -> list[list[str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    commands: list[list[str]] = []
    for source in sorted(gpd_dir.glob("*.gpd")):
        normalized = out_dir / f"{source.stem}.genePred"
        target = out_dir / f"{source.stem}.gtf"
        _normalize_idp_gpd(source, normalized)
        command = [str(gene_pred_to_gtf), "file", str(normalized), str(target)]
        commands.append(command)
        if not dry_run:
            subprocess.run(command, check=True)
    return commands


def run_suppa(samples: list[Sample], gtf_dir: Path, out_dir: Path, suppa: str, dry_run: bool) -> list[list[str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    commands: list[list[str]] = []
    for gtf in sorted(gtf_dir.glob("*.gtf")):
        prefix = out_dir / gtf.stem
        commands.append([suppa, "generateEvents", "-i", str(gtf), "-o", str(prefix), "-f", "ioe", "-e", "SE", "SS", "MX", "RI", "FL"])
        commands.append([suppa, "generateEvents", "-i", str(gtf), "-o", str(prefix), "-f", "ioi"])
        if samples and all(sample.tpm for sample in samples):
            # SUPPA expects a transcript expression matrix; for processed-input reruns this is commonly shared.
            tpm = str(samples[0].tpm)
            commands.append([suppa, "psiPerIsoform", "-g", str(gtf), "-e", tpm, "-o", str(prefix)])
    if not dry_run:
        for command in commands:
            subprocess.run(command, check=True)
    return commands


def summarize_outputs(results_dir: Path, out: Path) -> None:
    rows: list[dict[str, str | int]] = []
    for path in sorted(results_dir.rglob("*")):
        if path.suffix not in {".gpd", ".gtf", ".ioe", ".psi"}:
            continue
        rows.append(
            {
                "file": str(path),
                "type": path.suffix.lstrip("."),
                "records": _count_records(path),
            }
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "type", "records"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _format_idp_config(values: dict[str, Any]) -> str:
    ordered_keys = [
        "Nthread",
        "Njun_limit",
        "Niso_limit",
        "L_exon_limit",
        "L_min_intron",
        "Bfile_Npt",
        "Bfile_Nbin",
        "exon_construction_junction_span",
        "python_path",
        "temp_foldername",
        "output_foldername",
        "LR_gpd_pathfilename",
        "LR_psl_pathfilename",
        "LR_pathfilename",
        "genome_pathfilename",
        "aligner_choice",
        "gmap_index_pathfoldername",
        "blat_executable_pathfilename",
        "gmap_executable_pathfilename",
        "seqmap_executable_pathfilename",
        "ref_annotation_pathfilename",
        "allref_annotation_pathfilename",
        "SR_jun_pathfilename",
        "SR_sam_pathfilename",
        "read_length",
        "min_junction_overlap_len",
        "CAGE_data_filename",
        "detected_exp_len",
        "I_refjun_isoformconstruction",
        "I_ref5end_isoformconstruction",
        "I_ref3end_isoformconstruction",
        "estimator_choice",
        "FPR",
        "min_isoform_fraction",
        "min_isoform_rpkm",
    ]
    lines = ["# Generated by idp-beta make-idp-configs"]
    for key in ordered_keys:
        value = values.get(key, "")
        lines.append(f"{key} = {value}")
    return "\n".join(lines) + "\n"


def _long_read_config_key(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".gpd":
        return "LR_gpd_pathfilename"
    if suffix == ".psl":
        return "LR_psl_pathfilename"
    return "LR_pathfilename"


def _normalize_idp_gpd(source: Path, target: Path) -> None:
    with source.open() as inp, target.open("w") as out:
        for line in inp:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 11:
                fields = [fields[0], *fields[2:11]]
            out.write("\t".join(fields) + "\n")


def _count_records(path: Path) -> int:
    with path.open(errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def merge_expression_tables(idp_output_dir: Path, out: Path, samples: list[Sample], expression_name: str = "isoform.exp") -> None:
    matrix: dict[str, dict[str, str]] = {}
    sample_ids = [sample.sample_id for sample in samples]
    for sample in samples:
        values: dict[str, float] = {}
        for source in sorted(idp_output_dir.glob(f"{sample.sample_id}*/{expression_name}")):
            with source.open() as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.split()
                    if len(fields) < 2:
                        continue
                    values[fields[0]] = values.get(fields[0], 0.0) + float(fields[1])
        for transcript, value in values.items():
            matrix.setdefault(transcript, {})[sample.sample_id] = f"{value:g}"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", *sample_ids])
        for transcript in sorted(matrix):
            writer.writerow([transcript, *[matrix[transcript].get(sample_id, "0") for sample_id in sample_ids]])


def snakemake_command(
    samples: Path,
    references: Path,
    workflow: Path,
    snakefile: Path,
    cores: int,
    dry_run: bool,
    extra_args: list[str] | None = None,
) -> list[str]:
    command = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--cores",
        str(cores),
        "--config",
        f"samples={samples}",
        f"references={references}",
        f"workflow={workflow}",
    ]
    if dry_run:
        command.append("--dry-run")
    if extra_args:
        command.extend(extra_args)
    return command


def run_snakemake(command: list[str], dry_run: bool, allow_missing: bool = False) -> None:
    if shutil.which("snakemake") is None:
        message = "Snakemake is not installed. Command to run after installing Snakemake:\n" + " ".join(command)
        if dry_run or allow_missing:
            print(message)
            return
        raise ConfigError(message)
    subprocess.run(command, check=True)


def workflow_chromosomes(workflow: dict[str, Any]) -> list[str]:
    raw = workflow.get("chromosomes", DEFAULT_CHROMOSOMES)
    if isinstance(raw, str):
        return [item.strip() for item in raw.split(",") if item.strip()]
    if isinstance(raw, list):
        return [str(item) for item in raw]
    raise ConfigError("workflow.yml chromosomes must be a list or comma-separated string")


def _sample_long_read(sample: Sample, chromosome: str | None) -> Path:
    if sample.long_read:
        return sample.long_read
    if chromosome:
        return Path(f"results/long_read/gmap/{sample.sample_id}.{chromosome}.psl")
    return Path(f"results/fmlrc/{sample.sample_id}.corrected.fa")


def _sample_sr_sam(sample: Sample, chromosome: str | None) -> Path:
    if sample.sr_sam:
        return sample.sr_sam
    if chromosome:
        return Path(f"results/idp_inputs/{sample.sample_id}.{chromosome}.splicemap-like.sam")
    return Path(f"results/idp_inputs/{sample.sample_id}.splicemap-like.sam")


def _sample_junction_bed(sample: Sample, chromosome: str | None) -> Path:
    if sample.sr_junction_bed:
        return sample.sr_junction_bed
    if chromosome:
        return Path(f"results/idp_inputs/{sample.sample_id}.{chromosome}.splicemap-like.junctions.bed")
    return Path(f"results/idp_inputs/{sample.sample_id}.splicemap-like.junctions.bed")

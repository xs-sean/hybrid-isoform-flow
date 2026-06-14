from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

from .adapters import convert_sam_to_idp_inputs
from .config import ConfigError, read_samples, read_simple_yaml, validate_configs
from .fixture import write_project_templates, write_synthetic_fixture
from .pipeline import (
    gpd_to_gtf,
    make_idp_configs,
    merge_expression_tables,
    merge_gpd_outputs,
    merge_sample_gpd_outputs,
    run_idp_configs,
    run_snakemake,
    run_suppa,
    snakemake_command,
    summarize_outputs,
    workflow_chromosomes,
)


DEFAULT_SAMPLES = Path("config/samples.tsv")
DEFAULT_REFERENCES = Path("config/references.yml")
DEFAULT_IDP = Path("config/idp.yml")
DEFAULT_WORKFLOW = Path("config/workflow.yml")
DEFAULT_SNAKEFILE = Path("workflow/Snakefile")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="idp-beta")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("init", help="Write public workflow config templates.")
    p.add_argument("--out-dir", type=Path, default=Path("."))
    p.add_argument("--force", action="store_true")

    p = sub.add_parser("run", help="Run the Snakemake workflow.")
    add_workflow_args(p)
    p.add_argument("--cores", type=int, default=1)
    p.add_argument("snakemake_args", nargs=argparse.REMAINDER)

    p = sub.add_parser("dry-run", help="Print or run a Snakemake dry-run for the full workflow.")
    add_workflow_args(p)
    p.add_argument("--cores", type=int, default=1)
    p.add_argument("snakemake_args", nargs=argparse.REMAINDER)

    p = sub.add_parser("make-synthetic-fixture", help="Create a tiny synthetic fixture for tests and examples.")
    p.add_argument("--out-dir", type=Path, default=Path("examples/synthetic"))
    p.add_argument("--force", action="store_true")

    add_config_args(sub.add_parser("validate-config", help="Validate sample, reference, and IDP settings."))
    p = sub.add_parser("make-idp-configs", help="Generate one IDP run.cfg per sample.")
    add_config_args(p)
    p.add_argument("--out-dir", type=Path, default=Path("work/idp_configs"))
    p.add_argument("--sample-id", default="")
    p.add_argument("--chromosome", default="")

    p = sub.add_parser("run-idp", help="Run upstream IDP-fusion for generated configs.")
    p.add_argument("--config-dir", type=Path, default=Path("work/idp_configs"))
    p.add_argument("--config-file", type=Path)
    p.add_argument("--idp-root", type=Path, default=Path("third_party/IDP-fusion"))
    p.add_argument("--container-image", default="")
    p.add_argument("--workdir", type=Path, default=Path.cwd())
    p.add_argument("--dry-run", action="store_true")

    p = sub.add_parser("merge-gpd", help="Merge per-sample IDP GPD outputs.")
    add_config_args(p)
    p.add_argument("--workflow", type=Path, default=DEFAULT_WORKFLOW)
    p.add_argument("--idp-output-dir", type=Path, default=Path("results/idp"))
    p.add_argument("--out-dir", type=Path, default=Path("results/gpd"))
    p.add_argument("--per-sample", action="store_true")

    p = sub.add_parser("gpd-to-gtf", help="Convert merged GPD/genePred outputs to GTF.")
    add_config_args(p)
    p.add_argument("--gpd-dir", type=Path, default=Path("results/gpd"))
    p.add_argument("--out-dir", type=Path, default=Path("results/gtf"))
    p.add_argument("--dry-run", action="store_true")

    p = sub.add_parser("run-suppa", help="Generate SUPPA IOE/IOI/PSI commands from GTF and TPM files.")
    add_config_args(p)
    p.add_argument("--gtf-dir", type=Path, default=Path("results/gtf"))
    p.add_argument("--out-dir", type=Path, default=Path("results/suppa"))
    p.add_argument("--dry-run", action="store_true")

    p = sub.add_parser("summarize", help="Summarize GPD, IOE, and PSI outputs.")
    p.add_argument("--results-dir", type=Path, default=Path("results"))
    p.add_argument("--out", type=Path, default=Path("results/summary.tsv"))

    p = sub.add_parser("adapt-hisat2", help="Convert HISAT2 SAM/BAM to IDP-compatible SAM and junction BED.")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--out-sam", type=Path, required=True)
    p.add_argument("--out-bed", type=Path, required=True)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--chromosome", default="")
    p.add_argument("--min-support", type=int, default=1)

    p = sub.add_parser("merge-expression", help="Merge per-sample IDP isoform.exp files into a TPM-like matrix.")
    add_config_args(p)
    p.add_argument("--idp-output-dir", type=Path, default=Path("results/idp"))
    p.add_argument("--out", type=Path, default=Path("results/expression/allSamples.tpm.tsv"))

    args = parser.parse_args(argv)
    try:
        if args.command == "init":
            written = write_project_templates(args.out_dir, args.force)
            for path in written:
                print(path)
        elif args.command == "run":
            command = snakemake_command(args.samples, args.references, args.workflow, args.snakefile, args.cores, False, _clean_remainder(args.snakemake_args))
            run_snakemake(command, dry_run=False)
        elif args.command == "dry-run":
            validate_configs(args.samples, args.references, args.idp, check_files=False)
            command = snakemake_command(args.samples, args.references, args.workflow, args.snakefile, args.cores, True, _clean_remainder(args.snakemake_args))
            run_snakemake(command, dry_run=True, allow_missing=True)
        elif args.command == "make-synthetic-fixture":
            written = write_synthetic_fixture(args.out_dir, args.force)
            print(f"Wrote {len(written)} synthetic fixture file(s) to {args.out_dir}")
        elif args.command == "validate-config":
            for message in validate_configs(args.samples, args.references, args.idp, check_files=not args.no_check_files):
                print(message)
        elif args.command == "make-idp-configs":
            samples = read_samples(args.samples)
            references = read_simple_yaml(args.references)
            idp = read_simple_yaml(args.idp)
            written = make_idp_configs(samples, references, idp, args.out_dir, args.sample_id or None, args.chromosome or None)
            print(f"Wrote {len(written)} IDP config file(s) to {args.out_dir}")
        elif args.command == "run-idp":
            commands = run_idp_configs(args.config_dir, args.idp_root, args.container_image, args.workdir, args.dry_run, args.config_file)
            _print_commands(commands, args.dry_run)
        elif args.command == "merge-gpd":
            if args.per_sample:
                samples = read_samples(args.samples)
                workflow = read_simple_yaml(args.workflow)
                written = merge_sample_gpd_outputs(args.idp_output_dir, args.out_dir, samples, workflow_chromosomes(workflow))
            else:
                written = merge_gpd_outputs(args.idp_output_dir, args.out_dir)
            print(f"Wrote {len(written)} merged GPD file(s) to {args.out_dir}")
        elif args.command == "gpd-to-gtf":
            references = read_simple_yaml(args.references)
            commands = gpd_to_gtf(args.gpd_dir, args.out_dir, Path(str(references["genePredToGtf"])), args.dry_run)
            _print_commands(commands, args.dry_run)
        elif args.command == "run-suppa":
            samples = read_samples(args.samples)
            references = read_simple_yaml(args.references)
            commands = run_suppa(samples, args.gtf_dir, args.out_dir, str(references.get("suppa", "suppa.py")), args.dry_run)
            _print_commands(commands, args.dry_run)
        elif args.command == "summarize":
            summarize_outputs(args.results_dir, args.out)
            print(f"Wrote summary to {args.out}")
        elif args.command == "adapt-hisat2":
            stats = convert_sam_to_idp_inputs(args.input, args.out_sam, args.out_bed, args.read_length, args.chromosome or None, min_support=args.min_support)
            print(f"Kept {stats.kept_records}/{stats.total_records} SAM records and wrote {stats.junctions} junction(s).")
        elif args.command == "merge-expression":
            samples = read_samples(args.samples)
            merge_expression_tables(args.idp_output_dir, args.out, samples)
            print(f"Wrote merged expression matrix to {args.out}")
    except (ConfigError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


def add_config_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--samples", type=Path, default=DEFAULT_SAMPLES)
    parser.add_argument("--references", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--idp", type=Path, default=DEFAULT_IDP)
    parser.add_argument("--no-check-files", action="store_true", help="Validate schema only; do not require configured paths to exist.")


def add_workflow_args(parser: argparse.ArgumentParser) -> None:
    add_config_args(parser)
    parser.add_argument("--workflow", type=Path, default=DEFAULT_WORKFLOW)
    parser.add_argument("--snakefile", type=Path, default=DEFAULT_SNAKEFILE)


def _print_commands(commands: list[list[str]], dry_run: bool) -> None:
    if dry_run:
        for command in commands:
            print(" ".join(shlex.quote(part) for part in command))


def _clean_remainder(args: list[str]) -> list[str]:
    if args and args[0] == "--":
        return args[1:]
    return args


if __name__ == "__main__":
    raise SystemExit(main())

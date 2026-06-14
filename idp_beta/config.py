from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """Raised when an IDP_BETA configuration is invalid."""


@dataclass(frozen=True)
class Sample:
    sample_id: str
    condition: str
    isoseq_reads: Path | None = None
    short_read_1: Path | None = None
    short_read_2: Path | None = None
    external_tpm: Path | None = None
    long_read: Path | None = None
    sr_sam: Path | None = None
    sr_junction_bed: Path | None = None

    @property
    def tpm(self) -> Path | None:
        return self.external_tpm

    @property
    def has_raw_inputs(self) -> bool:
        return bool(self.isoseq_reads and self.short_read_1 and self.short_read_2)

    @property
    def has_processed_inputs(self) -> bool:
        return bool(self.long_read and self.sr_sam and self.sr_junction_bed)


def read_samples(path: Path) -> list[Sample]:
    if not path.exists():
        raise ConfigError(f"Missing sample sheet: {path}")
    with path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(_non_comment_lines(handle), delimiter="\t")]
    raw_required = {"sample_id", "condition", "isoseq_reads", "short_read_1", "short_read_2"}
    processed_required = {"sample_id", "condition", "long_read", "sr_sam", "sr_junction_bed"}
    if not rows:
        return []
    columns = set(rows[0])
    if not raw_required <= columns and not processed_required <= columns:
        raw_missing = raw_required - columns
        processed_missing = processed_required - columns
        raise ConfigError(
            f"{path} must use either raw workflow columns "
            f"({', '.join(sorted(raw_required))}) or processed-input columns "
            f"({', '.join(sorted(processed_required))}). Missing raw columns: "
            f"{', '.join(sorted(raw_missing)) or 'none'}; missing processed columns: "
            f"{', '.join(sorted(processed_missing)) or 'none'}"
        )
    seen: set[str] = set()
    samples: list[Sample] = []
    for row_num, row in enumerate(rows, start=2):
        sample_id = row["sample_id"].strip()
        condition = row["condition"].strip()
        if not sample_id:
            raise ConfigError(f"{path}:{row_num} has an empty sample_id")
        if sample_id in seen:
            raise ConfigError(f"{path}:{row_num} duplicates sample_id {sample_id!r}")
        if not condition:
            raise ConfigError(f"{path}:{row_num} has an empty condition")
        seen.add(sample_id)
        external_tpm_value = row.get("external_tpm", row.get("tpm", "")).strip()
        samples.append(
            Sample(
                sample_id=sample_id,
                condition=condition,
                isoseq_reads=_optional_path(row.get("isoseq_reads", "")),
                short_read_1=_optional_path(row.get("short_read_1", "")),
                short_read_2=_optional_path(row.get("short_read_2", "")),
                external_tpm=Path(external_tpm_value) if external_tpm_value else None,
                long_read=_optional_path(row.get("long_read", "")),
                sr_sam=_optional_path(row.get("sr_sam", "")),
                sr_junction_bed=_optional_path(row.get("sr_junction_bed", "")),
            )
        )
    return samples


def read_simple_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ConfigError(f"Missing config file: {path}")
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]
    prepared = _prepared_yaml_lines(path)
    for index, (line_num, line) in enumerate(prepared):
        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()
        if ":" not in stripped:
            raise ConfigError(f"{path}:{line_num} is not simple key/value YAML")
        key, value = stripped.split(":", 1)
        key = key.strip()
        value = value.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if value == "":
            next_indent = prepared[index + 1][1].find(prepared[index + 1][1].lstrip(" ")) if index + 1 < len(prepared) else -1
            if next_indent > indent:
                child: dict[str, Any] = {}
                parent[key] = child
                stack.append((indent, child))
            else:
                parent[key] = ""
        else:
            parent[key] = _parse_scalar(value)
    return root


def validate_configs(samples_path: Path, references_path: Path, idp_path: Path, check_files: bool = True) -> list[str]:
    messages: list[str] = []
    samples = read_samples(samples_path)
    references = read_simple_yaml(references_path)
    idp = read_simple_yaml(idp_path)
    if not samples:
        messages.append(f"No samples are defined in {samples_path}; add rows before running the workflow.")
    required_reference_keys = ["reference_gpd", "allref_gpd", "genome_fasta", "gmap_database", "hisat2_index", "chromosome_sizes", "genePredToGtf"]
    missing_refs = [key for key in required_reference_keys if not references.get(key)]
    if missing_refs:
        raise ConfigError(f"{references_path} is missing values for: {', '.join(missing_refs)}")
    if references.get("build") and str(references["build"]).lower() not in {"grch38", "hg38"}:
        raise ConfigError(f"{references_path} build must be GRCh38/hg38 for the public workflow")
    required_idp_keys = ["Nthread", "Njun_limit", "Niso_limit", "L_exon_limit", "L_min_intron", "read_length", "min_junction_overlap_len", "estimator_choice"]
    missing_idp = [key for key in required_idp_keys if key not in idp]
    if missing_idp:
        raise ConfigError(f"{idp_path} is missing values for: {', '.join(missing_idp)}")
    if check_files:
        paths: list[Path] = []
        for sample in samples:
            if sample.has_raw_inputs:
                paths.extend([sample.isoseq_reads, sample.short_read_1, sample.short_read_2])
            if sample.has_processed_inputs:
                paths.extend([sample.long_read, sample.sr_sam, sample.sr_junction_bed])
            if sample.external_tpm:
                paths.append(sample.external_tpm)
        paths.extend(Path(str(references[key])) for key in required_reference_keys)
        missing = [str(path) for path in paths if path is not None and not path.exists()]
        if missing:
            raise ConfigError("Missing configured files:\n" + "\n".join(f"  - {item}" for item in missing))
    messages.append(f"Validated {len(samples)} sample(s).")
    return messages


def _non_comment_lines(handle):
    for line in handle:
        if line.lstrip().startswith("#") or not line.strip():
            continue
        yield line


def _prepared_yaml_lines(path: Path) -> list[tuple[int, str]]:
    lines: list[tuple[int, str]] = []
    for line_num, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.split("#", 1)[0].rstrip()
        if line.strip():
            lines.append((line_num, line))
    return lines


def _parse_scalar(value: str) -> Any:
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if value in {"null", "None", "~"}:
        return None
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(item.strip()) for item in inner.split(",")]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _optional_path(value: str | None) -> Path | None:
    if value is None:
        return None
    stripped = value.strip()
    return Path(stripped) if stripped else None

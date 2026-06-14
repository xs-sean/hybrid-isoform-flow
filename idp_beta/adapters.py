from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


SUPPORTED_IDP_CIGAR = {"M", "N", "I", "D"}
CONSUMES_REFERENCE = {"M", "N", "D"}
CONSUMES_QUERY = {"M", "I"}


@dataclass(frozen=True)
class AdapterStats:
    total_records: int
    kept_records: int
    dropped_records: int
    junctions: int


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    if not cigar or cigar == "*":
        return []
    parts = re.findall(r"(\d+)([A-Z=])", cigar)
    if "".join(f"{length}{op}" for length, op in parts) != cigar:
        return []
    normalized: list[tuple[int, str]] = []
    for length, op in parts:
        op = "M" if op in {"=", "X"} else op
        normalized.append((int(length), op))
    return normalized


def cigar_to_idp(cigar: str) -> str | None:
    parts = parse_cigar(cigar)
    if not parts or any(op not in SUPPORTED_IDP_CIGAR for _, op in parts):
        return None
    return "".join(f"{length}{op}" for length, op in parts)


def query_length(parts: Iterable[tuple[int, str]]) -> int:
    return sum(length for length, op in parts if op in CONSUMES_QUERY)


def convert_sam_to_idp_inputs(
    sam_path: Path,
    out_sam: Path,
    out_bed: Path,
    read_length: int = 0,
    chromosome: str | None = None,
    flank: int = 50,
    min_support: int = 1,
) -> AdapterStats:
    out_sam.parent.mkdir(parents=True, exist_ok=True)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    kept = 0
    junction_reads: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    with _open_sam_text(sam_path) as inp, out_sam.open("w") as sam_out:
        for raw in inp:
            if raw.startswith("@") or not raw.strip():
                continue
            total += 1
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 11 or not _is_primary_mapped(fields):
                continue
            if chromosome and fields[2] != chromosome:
                continue
            idp_cigar = cigar_to_idp(fields[5])
            if idp_cigar is None:
                continue
            parts = parse_cigar(idp_cigar)
            if read_length and query_length(parts) != read_length:
                continue
            fields[5] = idp_cigar
            sam_out.write("\t".join(fields) + "\n")
            kept += 1
            for junction in _junctions_from_alignment(fields[2], int(fields[3]) - 1, parts):
                junction_reads[junction].add(fields[0])
    written_junctions = _write_junction_bed(junction_reads, out_bed, flank, min_support)
    return AdapterStats(total, kept, total - kept, written_junctions)


@contextmanager
def _open_sam_text(path: Path) -> Iterator[TextIO]:
    if path.suffix.lower() == ".bam":
        process = subprocess.Popen(["samtools", "view", "-h", str(path)], stdout=subprocess.PIPE, text=True)
        if process.stdout is None:
            raise RuntimeError("samtools did not provide stdout")
        try:
            yield process.stdout
        finally:
            process.stdout.close()
            returncode = process.wait()
            if returncode:
                raise subprocess.CalledProcessError(returncode, ["samtools", "view", "-h", str(path)])
    else:
        with path.open() as handle:
            yield handle


def _is_primary_mapped(fields: list[str]) -> bool:
    try:
        flag = int(fields[1])
    except ValueError:
        return False
    if fields[2] == "*" or fields[3] == "0":
        return False
    return not (flag & 0x4 or flag & 0x100 or flag & 0x800)


def _junctions_from_alignment(chromosome: str, start0: int, cigar_parts: list[tuple[int, str]]) -> list[tuple[str, int, int]]:
    ref_pos = start0
    junctions: list[tuple[str, int, int]] = []
    for length, op in cigar_parts:
        if op == "N":
            left = ref_pos
            ref_pos += length
            junctions.append((chromosome, left, ref_pos))
        elif op in CONSUMES_REFERENCE:
            ref_pos += length
    return junctions


def _write_junction_bed(
    junction_reads: dict[tuple[str, int, int], set[str]],
    out_bed: Path,
    flank: int,
    min_support: int,
) -> int:
    count = 0
    with out_bed.open("w") as handle:
        for (chromosome, left, right), reads in sorted(junction_reads.items()):
            support = len(reads)
            if support < min_support:
                continue
            start = max(0, left - flank)
            end = right + flank
            left_block = left - start
            right_block = end - right
            name = f"junction_[{support}]({support}/0)"
            row = [
                chromosome,
                str(start),
                str(end),
                name,
                "50",
                "+",
                str(start),
                str(end),
                "0,0,0",
                "2",
                f"{left_block},{right_block}",
                f"0,{right - left + left_block}",
            ]
            handle.write("\t".join(row) + "\n")
            count += 1
    return count

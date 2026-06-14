from __future__ import annotations

from pathlib import Path


PROJECT_TEMPLATES = {
    "config/samples.tsv": """sample_id\tcondition\tisoseq_reads\tshort_read_1\tshort_read_2\texternal_tpm
# Con1\tcontrol\tdata/isoseq/Con1.flnc.fa\tdata/illumina/Con1_R1.fastq.gz\tdata/illumina/Con1_R2.fastq.gz\tdata/expression/allSamples.tpm.tsv
""",
    "config/references.yml": """# Absolute paths are recommended for production reruns.
build: GRCh38
genome_fasta: data/references/GRCh38.fa
chromosome_sizes: data/references/GRCh38.chrom.sizes
hisat2_index: data/references/hisat2/GRCh38
gmap_database: data/references/gmap/GRCh38
reference_gpd: data/references/GRCh38.annotation.gpd
allref_gpd: data/references/GRCh38.all.gene.refFlat.txt
genePredToGtf: /usr/local/bin/genePredToGtf
suppa: suppa.py
blat_executable_pathfilename: blat
gmap_executable_pathfilename: gmap
seqmap_executable_pathfilename: seqmap
""",
    "config/workflow.yml": """chromosomes: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY, chrM]
threads: 10
read_length: 150
trimmomatic_jar: trimmomatic
trimmomatic_options: "TRAILING:20 MINLEN:235 CROP:235"
hisat2_options: "--dta"
gmap_options: "-f 1 -n 1"
# Available fmlrc_command placeholders: {input}/{long_read}, {output},
# {short_read_1}, {short_read_2}, {trimmed_read_1}, {trimmed_read_2},
# {short_read_msbwt}, and {fmlrc_options}.
fmlrc_command: "fmlrc {fmlrc_options} {short_read_msbwt} {input} {output}"
fmlrc_options: ""
short_read_msbwt: data/references/short_reads.msbwt.npy
container_image: idp-beta-py2
suppa_events: "SE SS MX RI FL"
""",
}


SYNTHETIC_FILES = {
    "config/samples.tsv": """sample_id\tcondition\tisoseq_reads\tshort_read_1\tshort_read_2\texternal_tpm
ToyCon\tcontrol\tdata/isoseq/ToyCon.flnc.fa\tdata/illumina/ToyCon_R1.fastq\tdata/illumina/ToyCon_R2.fastq\tdata/expression/allSamples.tpm.tsv
ToyIPF\tipf\tdata/isoseq/ToyIPF.flnc.fa\tdata/illumina/ToyIPF_R1.fastq\tdata/illumina/ToyIPF_R2.fastq\tdata/expression/allSamples.tpm.tsv
""",
    "config/references.yml": """build: GRCh38
genome_fasta: data/references/toy.GRCh38.fa
chromosome_sizes: data/references/toy.GRCh38.chrom.sizes
hisat2_index: data/references/hisat2/toy.GRCh38
gmap_database: data/references/gmap/toy.GRCh38
reference_gpd: data/references/toy.GRCh38.gpd
allref_gpd: data/references/toy.GRCh38.refFlat.txt
genePredToGtf: genePredToGtf
suppa: suppa.py
blat_executable_pathfilename: blat
gmap_executable_pathfilename: gmap
seqmap_executable_pathfilename: seqmap
""",
    "config/workflow.yml": """chromosomes: [chrToy]
threads: 1
read_length: 10
trimmomatic_jar: trimmomatic
trimmomatic_options: "MINLEN:10"
hisat2_options: "--dta"
gmap_options: "-f 1 -n 1"
fmlrc_command: "cp {input} {output}"
fmlrc_options: ""
short_read_msbwt: data/references/toy.msbwt.npy
container_image: idp-beta-py2
suppa_events: "SE SS MX RI FL"
""",
    "data/references/toy.GRCh38.fa": """>chrToy
AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCC
""",
    "data/references/toy.GRCh38.chrom.sizes": "chrToy\t60\n",
    "data/references/toy.GRCh38.gpd": "ToyGene\tToyTx\tchrToy\t+\t0\t60\t0\t60\t2\t0,40,\t20,60,\n",
    "data/references/toy.GRCh38.refFlat.txt": "ToyGene\tToyTx\tchrToy\t+\t0\t60\t0\t60\t2\t0,40,\t20,60,\n",
    "data/references/toy.msbwt.npy": "placeholder\n",
    "data/isoseq/ToyCon.flnc.fa": ">ToyCon_read1\nAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT\n",
    "data/isoseq/ToyIPF.flnc.fa": ">ToyIPF_read1\nAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT\n",
    "data/illumina/ToyCon_R1.fastq": "@ToyCon_r1\nAAAAAAAAAA\n+\nIIIIIIIIII\n",
    "data/illumina/ToyCon_R2.fastq": "@ToyCon_r2\nTTTTTTTTTT\n+\nIIIIIIIIII\n",
    "data/illumina/ToyIPF_R1.fastq": "@ToyIPF_r1\nAAAAAAAAAA\n+\nIIIIIIIIII\n",
    "data/illumina/ToyIPF_R2.fastq": "@ToyIPF_r2\nTTTTTTTTTT\n+\nIIIIIIIIII\n",
    "data/expression/allSamples.tpm.tsv": "transcript_id\tToyCon\tToyIPF\nToyTx\t10\t12\n",
    "README.md": """# Synthetic IDP_BETA Fixture

This tiny fixture is artificial and contains no patient data. It is intended for static validation and Snakemake dry-runs, not biological interpretation.
""",
}


def write_project_templates(out_dir: Path, force: bool = False) -> list[Path]:
    return _write_files(out_dir, PROJECT_TEMPLATES, force)


def write_synthetic_fixture(out_dir: Path, force: bool = False) -> list[Path]:
    return _write_files(out_dir, SYNTHETIC_FILES, force)


def _write_files(out_dir: Path, files: dict[str, str], force: bool) -> list[Path]:
    written: list[Path] = []
    for relative, content in files.items():
        path = out_dir / relative
        if path.exists() and not force:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        written.append(path)
    return written

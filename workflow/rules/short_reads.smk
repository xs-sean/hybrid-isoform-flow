rule fastqc:
    input:
        r1=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].short_read_1,
        r2=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].short_read_2
    output:
        "results/qc/fastqc/{sample}.done"
    threads: 1
    shell:
        "mkdir -p results/qc/fastqc && fastqc -o results/qc/fastqc {input.r1} {input.r2} && touch {output}"


rule trim_short_reads:
    input:
        r1=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].short_read_1,
        r2=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].short_read_2,
        qc="results/qc/fastqc/{sample}.done"
    output:
        r1="work/trimmed/{sample}_R1.fastq.gz",
        r2="work/trimmed/{sample}_R2.fastq.gz",
        unpaired1="work/trimmed/{sample}_R1.unpaired.fastq.gz",
        unpaired2="work/trimmed/{sample}_R2.unpaired.fastq.gz"
    threads: THREADS
    params:
        trimmer=lambda wildcards: WORKFLOW.get("trimmomatic_jar", "trimmomatic"),
        options=lambda wildcards: WORKFLOW.get("trimmomatic_options", "")
    shell:
        """
        mkdir -p work/trimmed
        if [ "{params.trimmer}" = "trimmomatic" ]; then
          trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.unpaired1} {output.r2} {output.unpaired2} {params.options}
        else
          java -jar {params.trimmer} PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.unpaired1} {output.r2} {output.unpaired2} {params.options}
        fi
        """


rule hisat2_align:
    input:
        r1="work/trimmed/{sample}_R1.fastq.gz",
        r2="work/trimmed/{sample}_R2.fastq.gz"
    output:
        "work/hisat2/{sample}.sam"
    threads: THREADS
    params:
        index=lambda wildcards: REFERENCES["hisat2_index"],
        options=lambda wildcards: WORKFLOW.get("hisat2_options", "--dta")
    shell:
        "mkdir -p work/hisat2 && hisat2 -p {threads} {params.options} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output}"


rule hisat2_to_idp_inputs:
    input:
        "work/hisat2/{sample}.sam"
    output:
        sam="results/idp_inputs/{sample}.{chrom}.splicemap-like.sam",
        bed="results/idp_inputs/{sample}.{chrom}.splicemap-like.junctions.bed"
    params:
        read_length=READ_LENGTH
    shell:
        "idp-beta adapt-hisat2 --input {input} --out-sam {output.sam} --out-bed {output.bed} --read-length {params.read_length} --chromosome {wildcards.chrom}"

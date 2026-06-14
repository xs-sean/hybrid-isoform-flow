rule fmlrc_correct_long_reads:
    input:
        lr=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].isoseq_reads,
        r1="work/trimmed/{sample}_R1.fastq.gz",
        r2="work/trimmed/{sample}_R2.fastq.gz"
    output:
        "results/fmlrc/{sample}.corrected.fa"
    params:
        cmd=lambda wildcards, input, output: str(WORKFLOW.get("fmlrc_command", "fmlrc {fmlrc_options} {short_read_msbwt} {input} {output}")).format(
            fmlrc_options=WORKFLOW.get("fmlrc_options", ""),
            short_read_msbwt=WORKFLOW.get("short_read_msbwt", ""),
            input=input.lr,
            long_read=input.lr,
            short_read_1=SAMPLE_BY_ID[wildcards.sample].short_read_1,
            short_read_2=SAMPLE_BY_ID[wildcards.sample].short_read_2,
            trimmed_read_1=input.r1,
            trimmed_read_2=input.r2,
            output=output[0],
        )
    shell:
        """
        mkdir -p results/fmlrc
        {params.cmd}
        """


rule gmap_align_long_reads:
    input:
        "results/fmlrc/{sample}.corrected.fa"
    output:
        "results/long_read/gmap/{sample}.{chrom}.psl"
    threads: THREADS
    params:
        db=lambda wildcards: Path(str(REFERENCES["gmap_database"])).name,
        db_dir=lambda wildcards: str(Path(str(REFERENCES["gmap_database"])).parent),
        options=lambda wildcards: WORKFLOW.get("gmap_options", "-f 1 -n 1")
    shell:
        "mkdir -p results/long_read/gmap && gmap -D {params.db_dir} -d {params.db} -t {threads} {params.options} {input} > {output}"

wildcard_constraints:
    kind="isoform_detection|isoform_prediction|isoform"


rule merge_gpd_by_sample:
    input:
        lambda wildcards: expand("results/idp/{sample}.{chrom}/{kind}.gpd", sample=wildcards.sample, chrom=CHROMOSOMES, kind=wildcards.kind)
    output:
        "results/gpd/{sample}.{kind}.gpd"
    shell:
        "mkdir -p results/gpd && cat {input} > {output}"


rule merge_expression:
    input:
        lambda wildcards: expand("results/idp/{sample}.{chrom}/isoform.exp", sample=SAMPLE_IDS, chrom=CHROMOSOMES)
    output:
        "results/expression/allSamples.tpm.tsv"
    shell:
        "idp-beta merge-expression --samples {SAMPLES_PATH} --idp-output-dir results/idp --out {output}"


rule gpd_to_gtf:
    input:
        "results/gpd/{sample}.{kind}.gpd"
    output:
        "results/gtf/{sample}.{kind}.gtf"
    params:
        gene_pred_to_gtf=lambda wildcards: REFERENCES["genePredToGtf"]
    shell:
        "idp-beta gpd-to-gtf --references {REFERENCES_PATH} --gpd-dir results/gpd --out-dir results/gtf"


rule suppa_events:
    input:
        "results/gtf/{sample}.{kind}.gtf"
    output:
        "results/suppa/{sample}.{kind}.events.done"
    params:
        suppa=lambda wildcards: REFERENCES.get("suppa", "suppa.py"),
        prefix=lambda wildcards: f"results/suppa/{wildcards.sample}.{wildcards.kind}",
        events=" ".join(SUPPA_EVENTS)
    shell:
        """
        mkdir -p results/suppa
        {params.suppa} generateEvents -i {input} -o {params.prefix} -f ioe -e {params.events}
        {params.suppa} generateEvents -i {input} -o {params.prefix} -f ioi
        touch {output}
        """


rule suppa_psi:
    input:
        gtf="results/gtf/{sample}.{kind}.gtf",
        events="results/suppa/{sample}.{kind}.events.done",
        expression=sample_external_tpm
    output:
        "results/suppa/{sample}.{kind}.psi"
    params:
        suppa=lambda wildcards: REFERENCES.get("suppa", "suppa.py"),
        prefix=lambda wildcards: f"results/suppa/{wildcards.sample}.{wildcards.kind}"
    shell:
        """
        {params.suppa} psiPerIsoform -g {input.gtf} -e {input.expression} -o {params.prefix}
        if [ -f {params.prefix}_isoform.psi ]; then
          mv {params.prefix}_isoform.psi {output}
        elif [ -f {params.prefix}.psi ]; then
          mv {params.prefix}.psi {output}
        else
          touch {output}
        fi
        """


rule summarize_results:
    input:
        expand("results/suppa/{sample}.isoform_detection.psi", sample=SAMPLE_IDS)
    output:
        "results/summary.tsv"
    shell:
        "idp-beta summarize --results-dir results --out {output}"

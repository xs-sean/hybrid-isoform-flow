rule make_idp_config:
    input:
        lr="results/long_read/gmap/{sample}.{chrom}.psl",
        sam="results/idp_inputs/{sample}.{chrom}.splicemap-like.sam",
        bed="results/idp_inputs/{sample}.{chrom}.splicemap-like.junctions.bed"
    output:
        "work/idp_configs/{sample}.{chrom}.run.cfg"
    shell:
        "idp-beta make-idp-configs --samples {SAMPLES_PATH} --references {REFERENCES_PATH} --idp {IDP_PATH} --sample-id {wildcards.sample} --chromosome {wildcards.chrom} --out-dir work/idp_configs"


rule run_idp:
    input:
        cfg="work/idp_configs/{sample}.{chrom}.run.cfg"
    output:
        detection="results/idp/{sample}.{chrom}/isoform_detection.gpd",
        prediction="results/idp/{sample}.{chrom}/isoform_prediction.gpd",
        isoform="results/idp/{sample}.{chrom}/isoform.gpd",
        expression="results/idp/{sample}.{chrom}/isoform.exp"
    params:
        container=lambda wildcards: WORKFLOW.get("container_image", "")
    shell:
        """
        mkdir -p results/idp work/idp
        if [ -n "{params.container}" ]; then
          idp-beta run-idp --config-file {input.cfg} --container-image {params.container}
        else
          idp-beta run-idp --config-file {input.cfg}
        fi
        """

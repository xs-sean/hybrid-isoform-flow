# Workflow

The Snakemake workflow in `workflow/Snakefile` runs the public IDP_BETA pipeline from IsoSeq-derived long-read files and paired Illumina FASTQs to SUPPA PSI outputs.

Run from the repository root:

```bash
idp-beta dry-run
idp-beta run --cores 10
```

The workflow is split into rule modules:

- `rules/short_reads.smk`: FastQC, trimming, HISAT2, and IDP short-read adapters
- `rules/long_reads.smk`: FMLRC correction and GMAP long-read alignment
- `rules/idp.smk`: IDP config generation and IDP-fusion execution
- `rules/postprocess.smk`: GPD merging, GTF conversion, SUPPA, PSI, and summaries

The default workflow uses GRCh38 consistently. If a legacy analysis used mixed reference names or builds, normalize all inputs and annotations before running this public workflow.

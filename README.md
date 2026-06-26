# nf-core-umidedup

## Introduction

This repos is serving the development of subworkflow -- vidrl_umi_dedup
you could find this subworkflow in `subworkflows/local/vidrl_umi_dedup`

## Implement the subworkflow in your nextflow/nf-core pipeline

```
    # copy the vidrl_umi_dedup folder into subworkflows/local folder in your pipeline
    # nf-core modules install all the dependency in the subworkflows
    include { VIDRL_UMI_DEDUP } from './subworkflows/local/vidrl_umi_dedup/main'

    #the example for input in the main.nf in the repos
```

## Test Run

```
    git clone https://github.com/abcdtree/nf-core-umidedup.git
    cd nf-core-umidedup
    nextflow run nf-core-umidedup/main.nf --input "nf-core-umidedup/tests/data/fastq/*_R{1,2}.fastq.gz" --fasta nf-core-umidedup/tests/data/ref/ref.fasta -profile docker --outdir tests
```
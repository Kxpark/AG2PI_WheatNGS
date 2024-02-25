###################################################
#           Index for BWA-MEM2
###################################################

rule bwa_Demo_index:
    input: 
        "/data/ref/{Genome}.gz"
    output:
        "data/ref/{Genome}.0123",
        "data/ref/{Genome}.amb",
        "data/ref/{Genome}.ann",
        "data/ref/{Genome}.bwt.2bit.64",
        "data/ref/{Genome}.pac"
    conda:
        "../envs/IndexReference.yml"
    threads: config['BWA-index']['threads']
    params:
        prefix="data/ref/{Genome}"

    shell:
        "bwa-mem2 index -p {params.prefix} {input}"
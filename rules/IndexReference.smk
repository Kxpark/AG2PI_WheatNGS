###################################################
#  Download Reference Genome and Index for BWA-MEM2
###################################################
rule Download_Reference:
    input:
        config['Reference_Genome']['Url']
    output:
        "/data/ref/{file}.zip"
    conda:
        "../envs/IndexReference.yml"
    shell:
        "wget {input}"

rule unzip:
    input:
        "/data/ref/{file}.zip"
    output:
        "/data/ref/{file}.fa"
    conda:
        "../envs/IndexReference.yml"
    shell:
        "unzip {input}"


#56 core and 420 GB?
rule bwa_index:
    input: 
        "{Genome}.fa"
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

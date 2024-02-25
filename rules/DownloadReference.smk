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

rule Gzip:
    input:
        "/data/ref/{file}.fa"
    output:
        "/data/ref/{file}.fa.gz"
    conda:
        "../envs/IndexReference.yml"
    shell:
        "gzip {input}"
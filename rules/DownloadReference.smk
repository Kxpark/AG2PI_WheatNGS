###################################################
#  Download Reference Genome and Index for BWA-MEM2
###################################################
rule Download_Reference:
    output:
        expand("data/ref/{file}",file=config['Reference_Genome']['Zipfile'])
    conda:
        "../envs/IndexReference.yml"
    params:
        url=config['Reference_Genome']['Url']
    shell:
        "wget -P data/ref {params.url}"

rule unzip:
    input:
        expand("data/ref/{file}",file=config['Reference_Genome']['Zipfile'])
    output:
        config['Reference_Genome']['ref']
    conda:
        "../envs/IndexReference.yml"
    shell:
        "unzip -d data/ref {input}"
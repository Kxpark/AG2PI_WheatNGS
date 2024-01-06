##############################################
#  Pull Fastq files from SRA number
#############################################

rule SRAprefetch:
    input: 
        samples=config['SRA']
    output:
        "/data/SRA/{samples}/{samples}.sra"
    conda:
        "../envs/SRApull.yml"
    params:
        maxsize=config['SRA-ToolKit']['max-size']
    shell:
        "prefetch --output-directory /data/SRA {input} --max-size {maxsize}"

rule SRApull:
    input:
        "/data/SRA/{samples}/{samples}.sra"
    output:
        "/data/samples/{samples}_1.fastq",
        "/data/samples/{samples}_2.fastq"
    conda:
        "../envs/SRApull.yml"
    threads: config['SRA-ToolKit']['threads']
    shell:
        "fasterq-dump {input} -O /data/samples --split-files -e {threads}"

rule pzip:
    input:
        "/data/samples/{samples}{strand}.fastq"
    output:
        "/data/samples/{samples}{strand}.fastq.gz"
    conda:
        "../envs/SRApull.yml"
    threads: config['SRA-ToolKit']['threads']
    shell:
        "pigz -p {threads} {input}"

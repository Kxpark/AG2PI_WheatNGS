##############################################
#  Pull Fastq files from SRA number
#############################################

rule SRAtotal:
    output:
        "data/samples/{sample}_1.fastq",
        "data/samples/{sample}_2.fastq"
    conda:
        "../envs/SRApull.yml"
    params:
        maxsize=config['SRA-ToolKit']['max-size']
    threads: config['SRA-ToolKit']['threads']
    wildcard_constraints:
        sample="|".join(config['SRA']['samples'])
    shell:
        """
        prefetch {wildcards.sample} --output-directory data/SRA/ --max-size {params.maxsize} && \
        fasterq-dump data/SRA/{wildcards.sample}/{wildcards.sample}.sra -O ./data/samples --split-files -e {threads} 
        """

rule pzip:
    input:
        "data/samples/{sample}_{num}.fastq"
    output:
        "data/samples/{sample}_{num}.fastq.gz"
    conda:
        "../envs/SRApull.yml"
    threads: config['SRA-ToolKit']['threads']
    params:
        name=config['SRA']['samples']
    shell:
        "pigz -p {threads} {input}"

rule Done_SRA:
    input:
        expand("data/samples/{sample}_{num}.fastq.gz",sample=config['SRA']['samples'],num=[1,2])
    output:
        "/data/samples/Done.txt"
    shell:
       "touch {output}"
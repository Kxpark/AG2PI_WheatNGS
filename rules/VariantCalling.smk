##################################################
# BWA-mem2 to BCFtools Variant Calling Pipeline
##################################################
rule Align_reads:
    input:
        fastq1="data/samples/{sample}_1.fastq.gz",
        fastq2="data/samples/{sample}_2.fastq.gz"
    output:
        "Aligned_Sorted_BAM/{sample}.bam"
    conda:
        "../envs/VariantCalling.yml"    
    log:
        "logs/Align_reads/{sample}.log"
        
    threads: 11
    params:
        index=config["Reference_Genome"]['index']
    shell:
        "(bwa-mem2 mem -t {threads} {params.index} {input.fastq1} {input.fastq2}  | samtools view -b | samtools sort -o {output}) 2> {log}"

rule Remove_multiple_Alignment: 
    input:
        "Aligned_Sorted_BAM/{sample}.bam"
    output:
        temp("Single_Read_BAM/{sample}.bam")
    conda:
        "../envs/VariantCalling.yml" 
    shell:
        "samtools view -h {input} | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b -o {output}"

rule Mark_Dupplicates:
    input:
        "Single_Read_BAM/{sample}.bam"
    output:
        temp("Marked_Dup/{sample}.bam")
    conda:
        "../envs/VariantCalling.yml" 
    log:
        "logs/Dupstats/Dupstats.{sample}.txt"
    shell:
        "picard MarkDuplicates I={input} O={output} M={log}"

rule Read_Groups:
    input:
       "Marked_Dup/{sample}.bam"
    output:
        temp("ReadGroup_Bam/{sample}.bam")
    conda:
        "../envs/VariantCalling.yml" 
    params:
        name=config["Population"]
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} RGLB={params.name} RGPL=ILLUMINA RGPU=unit RGSM={wildcards.sample}"

rule Unmapped_removal:
    input:
        "ReadGroup_Bam/{sample}.bam"
    output:
        "Done_Bam/{sample}.bam"
    conda:
        "../envs/VariantCalling.yml" 
    shell:
        "samtools view -h -F4 {input} -o {output}" 

rule samtools_index:
    input:
        "Done_Bam/{sample}.bam"
    output:
        "Done_Bam/{sample}.bam.csi"
    conda:
        "../envs/VariantCalling.yml"  
    shell:
        "samtools index -c {input}"

rule coverage:
    input:
        "Done_Bam/{sample}.bam"
    output:
        "Coverage/{sample}.cov"
    conda:
        "../envs/VariantCalling.yml" 
    shell:
        "samtools coverage {input} > {output}"


rule Split_Chromosome:
    input:
        bam="Done_Bam/{sample}.bam",
        index="Done_Bam/{sample}.bam.csi"
    output:
        "Split_Bam/{sample}.{Chr}.bam"
    params:
        Chr=config['Chromosomes']
    conda:
        "../envs/VariantCalling.yml" 
    shell:
        "samtools view -b -o {output} {input.bam} {wildcards.Chr}"

def sraSamples2(wildcards):
    if config['Workflow_Settings']['SRApull']:
        return ["Split_Bam/"+file + ".{Chr}.bam" for file in config['SRA']['samples']]
    else:
        return ["Split_Bam/"+file + ".{Chr}.bam" for file in config['samples']]

rule inputgenerator:
    input: 
        sraSamples2
    output:
        "Bamlists/Bamlist.{Chr}.txt"
    run:
        for f in input:
            with open(output[0],'a') as out:
                out.write(f + '\n') 

rule Mpileup:
    input:
        "Bamlists/Bamlist.{Chr}.txt"            
    output:
        "Split_VCF/{Chr}.raw.vcf.gz"
    params:
        ref=config["Reference_Genome"]['ref']
    conda:
        "../envs/VariantCalling.yml" 
    shell:
        "bcftools mpileup -b {input} -f {params.ref} -q 20 -Ou | bcftools call -Oz -mv -o {output}"


rule Concat_Mpileup:
    input:
        expand("Split_VCF/{Chr}.raw.vcf.gz",Chr=config['Chromosomes'])

    output:
        expand("Final_VCF/{name}.raw.vcf.gz",name=config["Population"])
    conda:
        "../envs/VariantCalling.yml" 
    threads: 48
    shell:
        "bcftools concat --threads {threads} -Oz -o {output} {input}"
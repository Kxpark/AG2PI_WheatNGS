################################
#  Filtering/Imputing VCF #not homozygous
################################
def get_vcfs(wildcards):
   return expand("Filt_VCF/{name}.missing{fmissing}_MAF{maf}_imputed.vcf.gz",name=config["Population"],fmissing=config['VCF_Filtering']['Percent_Missing']
                ,maf=config['VCF_Filtering']['Minor_allele_frequency'])

rule VCF_filt_1:
    input:
        "Final_VCF/{name}.raw.vcf.gz"
    output:
        "Final_VCF/Filter_VCF/{name}.missing{fmissing}.vcf.gz"
    conda:
        "../envs/VCFprocessing.yml" 
    params:
        name=config["Population"],
        fmissing=config['VCF_Filtering']['Percent_Missing'],
        Qual=config['VCF_Filtering']['Minimum_Quality']
    shell:
        "bcftools view -i 'F_MISSING<{params.fmissing} && QUAL>={params.Qual}' -Oz -o {output} {input}"

rule Beagle_Impute:
    input:
        "Final_VCF/Filter_VCF/{name}.missing{fmissing}.vcf.gz"
    output:
        "Final_VCF/Imputed/{name}.missing{fmissing}_imputed.vcf.gz"
    conda:
        "../envs/VCFprocessing.yml" 
    params:
        prefix='Final_VCF/Imputed/{name}.missing{fmissing}_imputed'    
    shell:
        "beagle gt={input} out={params.prefix}"        

rule VCF_filt_2:
    input:
        "Final_VCF/Imputed/{name}.missing{fmissing}_imputed.vcf.gz"
    output:
        "Filt_VCF/{name}.missing{fmissing}_MAF{maf}_imputed.vcf.gz"
    conda:
        "../envs/VCFprocessing.yml" 
    params:
        maf=config['VCF_Filtering']['Minor_allele_frequency']
    log:
        "logs/VCFfiltering/Filt2{name}.missing{fmissing}_{maf}.log"
    shell:
        "bcftools view -i 'MAF[0]>{params.maf}' -Oz -o {output} {input} 2> {log}"

rule Done:
    input:
        get_vcfs
    output:
        "Done.txt"
    script:
       "scripts/Done.py"
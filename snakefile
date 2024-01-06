from snakemake.utils import min_version
min_version('7.3.0')

configfile: "config.yaml"

##############################################
#  Load Rules
##############################################
#include:"rules/Background.smk"

if config['Workflow_Settings']['IndexReference']:
    include: "rules/IndexReference.smk"
if config['Workflow_Settings']['SRApull']:
    include: "rules/SRApull.smk"
include: "rules/VariantCalling.smk"

if config['Workflow_Settings']['VariantFiltering']:
    
    rule all:
        input:
            "Done.txt"
    
    if config['Workflow_Settings']['Homozygous']:
            include: "rules/VCFprocessingHomozygous.smk"
    else:
            include:"rules/VCFprocessing.smk"
else:
    rule all:
        input:
            expand("Final_VCF/{name}.raw.vcf.gz",name=config["Population"])

# PreWorkflow Checkup!
onstart:
    try:
        needed_files=[]
        print("Checking to make sure all required files are here....")
        if not config['Workflow_Settings']['IndexReference']:
            needed_files= needed_files + ["data/ref/*.0123",
                            "data/ref/*.amb",
                            "data/ref/*ann",
                            "data/ref/*.bwt.2bit.64",
                            "data/ref/*.pac"]
        if not config['Workflow_Settings']['SRApull']:
            for sample in config['samples']:
                file= sample + '.fastq.gz'
                needed_files= needed_files + [file] 
        for filename in needed_files:
            if not os.path.exists(filename):
                missing=filename
    except:
        print("This file was not found "+missing )
        sys.exit(1)
    else:
        print('/n All files are there, lets Go!')

#Success or Error messages

onsuccess:
	print("Success! The Variant Calling workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")

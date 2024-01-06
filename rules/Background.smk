#
# Start-up General 
#
# =================================================================================================
#     Import Libraries
# =================================================================================================

#import glob
#import os
#import json
##import pandas as pd
#import socket, platform
#from snakemake.io import expand
#from snakemake.utils import R
#from snakemake.utils import min_version
#from snakemake.utils import validate
#from snakemake.io import glob_wildcards
#import re
#rom os.path import join, basename, dirname
#import pathlib
#rom os import path
#from datetime import datetime


def get_vcfs(wildcards):
   return expand("Filt_VCF/{name}.missing{fmissing}_MAF{maf}_imputed.vcf.gz",name=config["Population"],fmissing=config['VCF_Filtering']['Percent_Missing']
                ,maf=config['VCF_Filtering']['Minor_allele_frequency'])


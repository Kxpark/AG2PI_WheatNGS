# Wheat NGS Variant Calling Snakemake Workflow <img align = "right" width = "220" src="https://github.com/Kxpark/AG2PI_WheatNGS/assets/104218928/8d2578cd-ed79-480c-842b-5019fb29a725">

A reproducible snakemake workflow that can be customized to fit your Next Gen-Sequencing projects.

First things first, This work is supported by funding from the Agricultural Genome to Phenome Initiative (AG2PI),
 which is funded by USDA-NIFA award 2022-70412-3845. <img align = "bottom" width = "110" src="https://github.com/Kxpark/AG2PI_WheatNGS/assets/104218928/1a214335-e1f1-407f-8d13-9b27ffb848a0">

This the first draft of the workflow, so feedback is welcome! Over time it will just get better.

## Table of Contents

* [Summary](#summary)
* [Installation](#installation)


## Summary

**Howdy!** This is a [Snakemake](https://snakemake.github.io) pipeline developed to aid in reproducible variant discovery. Alignment of NGS data to a reference genome can be meticulous and when performing analysis on large populations or working with high-coverage whole genome sequencing running analysis in parallel is a must. However, it can be difficult to ensure that all of your scripts are standardized and formatted for each sample. This workflow will take the difficulties of ensuring all of the intermediate jobs run correctly. While this was designed with hexaploid wheat in mind, any species with a reference genome and sequence data can be used. 

This workflow can be used with personal NGS sequence data or scripts are included to use SRA-toolkit to download sequence data from NCBI. Options in the Config file can be adjusted to best fit your project goals.

**Current issues** 
  1) With the size of the wheat genome, I am currently having issues uploading sample data or a reference genome that is already indexed (as the .tar file is over 80GB)
  2) Along this issue even a very low-coverage GBS example sample line is too big. These issues will hopefully be resolved soon.
  3) A small issue with the SRA pulling script will hopefully be resolved very soon.


## Installation

**1. Clone the git repository** 
```{bash}
git clone https://github.com/Kxpark/AG2PI_WheatNGS.git
```

**Do you already have Anaconda or Bioconda or Miniforge on your system?**
  If so skip to the next step.
**Otherwise** if you are on a fresh Linux environment, it is recommended to download Miniforge first for best results.
```{bash}
wget https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Miniforge3-23.3.1-1-Linux-x86_64.sh
bash Miniforge3-23.3.1-1-Linux-x86_64.sh
```
Then follow the instructions for installation.

**2. Setting up the Environment**
Once you have miniforge/conda downloaded go ahead and build your Snakemake Environment.
```{bash }
cd AG2PI_WheatNGS/
conda env create -f environment.yml -n VariantCallingSnakemake
conda activate VariantCallingSnakemake
```
If you have issues you can create a new conda environment, and make sure to load Snakemake, and Mamba packages. 
Then the nested environments inside of the /envs directory will be autoloaded.

**3. Set up the Config file**

Make sure to change the settings in the Config file to match what you are looking for and to match the data you are using.

**4. Dry Run** 
Once you think you have everything set, go ahead and do a dry run. This lets Snakemake do a check all files and code are ready before actually running jobs.
```{bash }
snakemake -n --use-conda
```

**4. Run the Pipeline**

```{bash }
snakemake --cores 11 --use-conda
```

Adding the `--use-conda` flag allows snakemake to build the individual dependencies for each rule. Likewise, you can adjust the `--cores` flag to the number of cores you are working with,
or use --cores all.

if you want to specify a config file location use the `--configfile` flag to go to another `config.yaml` file.
```{bash }
snakemake --cores 11 --use-conda --configfile <path/to/config.yaml>
```


If you have other questions the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) is very helpful.

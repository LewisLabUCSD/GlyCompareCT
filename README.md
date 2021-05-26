![Logo](logo_name.png)
# GlyCompare_CT
GlyCompare command-line tool。 GlyCompare is a novel method wherein glycans from glycomic data are decomposed to a minimal set of intermediate substructures, thus incorporating shared intermediate glycan substructures into all comparisons of glycans. 

## Citation

Bokan Bao+, Benjamin P. Kellman+, Austin W. T. Chiang, Austin K. York, Mahmoud A. Mohammad, Morey W. Haymond, Lars Bode, and Nathan E. Lewis. 2019. “**Correcting for Sparsity and Non-Independence in Glycomic Data through a System Biology Framework.**” bioRxiv. https://doi.org/10.1101/693507.

## Installation

First, please make sure you have `conda` installed. Version recommendation: conda 4.9.2 and later versions. 
- Install `conda` on Windows: https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html
- Install `conda` on Mac OS: https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html

Please ```git clone``` the main branch to your target local directory. 
```bash
# get the repo
git clone https://github.com/yuz682/GlyCompare_CT.git
# enter the repo
cd GlyCompare_CT
```

All dependencies required to run GlyCompare_CT can be installed using `environment.yml`. A new conda environment is created with all dependencies installed. This step will take a while (10 - 15 minutes). 
```bash
# Create the environment with all required dependencies installed.
conda env create -f environment.yml
```

Activate the new environment `glyCompare_CT_env`. Then the preprocessing is all done.
```bash
# Activate conda environment
conda activate glyCompare_CT_env
```

## User manual


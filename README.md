[![Build Status](https://travis-ci.org/cytham/nanovar.svg?branch=master)](https://travis-ci.org/cytham/nanovar)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/nanovar)](https://pypi.org/project/nanovar/)
[![PyPI versions](https://img.shields.io/pypi/v/nanovar)](https://pypi.org/project/nanovar/)
[![Conda](https://img.shields.io/conda/v/bioconda/nanovar)]()
[![Github release](https://img.shields.io/github/v/release/cytham/nanovar?include_prereleases)](https://github.com/cytham/nanovar/releases)
[![PyPI license](https://img.shields.io/pypi/l/nanovar)](https://github.com/cytham/nanovar/blob/master/LICENSE.txt)

<p align="center">
  <img src="http://benoukraf-lab.com/wp-content/uploads/2019/05/Nanovarlogo.png" width="200" alt="accessibility text" align='left'>
</p>

# NanoVar - Structural variant caller using low-depth long-read sequencing

NanoVar is a neural-network-based genomic structural variant (SV) caller that utilizes low-depth long-read sequencing such as
 Oxford Nanopore Technologies (ONT). It characterizes SVs with high accuracy and speed using only 4x depth
  sequencing for homozygous SVs and 8x depth for heterozygous SVs. NanoVar reduces sequencing cost and computational requirements
   which makes it compatible with large cohort SV-association studies or routine clinical SV investigations.  


### Basic capabilities
* Performs long-read mapping (HS-Blastn, Chen et al., 2015) and SV discovery in a single rapid pipeline.
* Accurately characterizes SVs using long sequencing reads (High SV recall and precision in simulation datasets, overall F1
 score >0.9)  
* Characterizes six classes of SVs including novel-sequence insertion, deletion, inversion, tandem duplication, sequence
 transposition and translocation.  
* Requires 4x and 8x sequencing depth for detecting homozygous and heterozygous SVs respectively.  
* Rapid computational speed (Takes <3 hours to map and analyze 12 gigabases datasets (4x) using 24 CPU threads)  
* Approximates SV genotype


## Getting Started

### Operating system: 
* Linux (x86_64 architecture, tested in Ubuntu 14.04, 16.04, 18.04)  

### Installation:
There are three ways to install NanoVar:
#### Option 1: Conda (Recommended)
```
# Installing from bioconda automatically installs all dependencies 
conda install -c bioconda nanovar
```
#### Option 2: Pip (See dependencies below)
```
# Installing from PyPI requires own installation of dependencies, see below
pip3 install nanovar
```
#### Option 3: GitHub (See dependencies below)
```
# Installing from GitHub requires own installation of dependencies, see below
git clone https://github.com/cytham/nanovar.git 
cd nanovar 
pip install .
```
### Installation of dependencies
* bedtools >=2.26.0
* makeblastdb and windowmasker
* hs-blastn

Please make sure each executable binary is in PATH.
##### 1. _bedtools_
Please visit [here](https://bedtools.readthedocs.io/en/latest/content/installation.html) for instructions to install.

##### 2. _makeblastdb_ and _windowmasker_
```
# Download NCBI-BLAST v2.3.0+ from NCBI FTP server
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz

# Extract tar.gz
tar zxf ncbi-blast-2.3.0+-x64-linux.tar.gz

# Copy makeblastdb and windowmasker binaries to PATH (e.g. ~/bin)
cp ncbi-blast-2.3.0+/bin/makeblastdb ~/bin && cp ncbi-blast-2.3.0+/bin/windowmasker ~/bin
```
##### 2. _hs-blastn_
```
# Download and compile
git clone https://github.com/chenying2016/queries.git
cd queries/hs-blastn-src/
make

# Copy hs-blastn binary to path (e.g. ~/bin)
cp hs-blastn ~/bin
```
### Quick run

```
nanovar [Options] -t 24 -f hg38 read.fa ref.fa working_dir 
```

| Parameter | Argument | Comment |
| :--- | :--- | :--- |
| `-t` | num_threads | Indicate number of CPU threads to use |
| `-f` | gap_file | Choose built-in gap BED file to exclude gap regions in the reference genome. Built-in gap files include: hg19, hg38 and mm10 (Optional)|
| - | read.fa | Input long-read FASTA/FASTQ file |
| - | ref.fa | Input reference genome in FASTA format |
| - | working_dir | Specify working directory |


## Documentation
See [Wiki](https://github.com/cytham/nanovar/wiki) for more information.

## Versioning
See [CHANGELOG](https://github.com/cytham/nanovar/blob/master/CHANGELOG.txt)

## Citation
NanoVar: Accurate Characterization of Patients’ Genomic Structural Variants Using Low-Depth Nanopore Sequencing (Tham. et al, 2019)
https://www.biorxiv.org/content/10.1101/662940v1
## Authors

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)
* **Roberto Tirado Magallanes** - [rtmag](https://github.com/rtmag)
* **Touati Benoukraf** - [benoukraflab](https://github.com/benoukraflab)

## License

This project is licensed under GNU General Public License - see [LICENSE.txt](https://github.com/cytham/nanovar/blob/master/LICENSE.txt) for details.

## Simulation datasets
SV-simulated datasets used for evaluating SV calling accuracy can be downloaded [here](https://doi.org/10.5281/zenodo.2599376).

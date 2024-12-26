## Please note: Current v1.8.1 not compatible with Tensorflow >= 2.16.0, please downgrade to 2.15.1

`pip install tensorflow-cpu==2.15.1`

Please see issue [here](https://github.com/cytham/nanovar/issues/77)

We are actively working on this, thank you for your understanding.

## NanoVar - Structural variant caller using low-depth long-read sequencing
[![Build Status](https://app.travis-ci.com/cytham/nanovar.svg?branch=master)](https://app.travis-ci.com/github/cytham/nanovar)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/nanovar)](https://pypi.org/project/nanovar/)
[![PyPI versions](https://img.shields.io/pypi/v/nanovar)](https://pypi.org/project/nanovar/)
<!--[![Conda](https://img.shields.io/conda/v/bioconda/nanovar)](https://anaconda.org/bioconda/nanovar)-->
[![Github release](https://img.shields.io/github/v/release/cytham/nanovar?include_prereleases)](../../releases)
[![PyPI license](https://img.shields.io/pypi/l/nanovar)](./LICENSE.txt)
  
NanoVar is a genomic structural variant (SV) caller that utilizes low-depth long-read sequencing such as
 Oxford Nanopore Technologies (ONT). It characterizes SVs with using only 4x depth
  sequencing for homozygous SVs and 8x depth for heterozygous SVs.

### Basic capabilities
* Performs long-read mapping (Minimap2) and SV discovery in a single pipeline.
* Accurately characterizes SVs using long sequencing reads (High SV recall and precision in simulation datasets, overall F1
 score >0.9)  
* Characterizes six classes of SVs including novel-sequence insertion, deletion, inversion, tandem duplication, sequence
 transposition (TPO) and translocation (TRA).  
* Requires 4x and 8x sequencing depth for detecting homozygous and heterozygous SVs respectively.  
* Rapid computational speed (Takes <3 hours to map and analyze 12 gigabases datasets (4x) using 24 CPU threads)  
* Approximates SV genotype
* Identifies full-length LINE and SINE insertions (Marked by "TE=" in the INFO column of VCF file)
* Repeat element INS annotation using [NanoINSight](https://github.com/AsmaaSamyMohamedMahmoud/NanoINSight)
<!--
* Detect large chromosomal copy-number variation using [CytoCAD](https://github.com/cytham/cytocad)
| `--cnv` | hg38 | Perform large CNV detection using CytoCAD (Only works for hg38 genome)
-->

## Getting Started

### Quick run

```
nanovar [Options] -t 24 -f hg38 sample.fq/sample.bam ref.fa working_dir 
```

| Parameter | Argument | Comment |
| :--- | :--- | :--- |
| `-t` | num_threads | Indicate number of CPU threads to use |
| `-f` (Optional) | gap_file (Optional) | Choose built-in gap BED file or specify own file to exclude gap regions in the reference genome. Built-in gap files include: hg19, hg38 and mm10 |
| - | sample.fq/sample.bam | Input long-read FASTA/FASTQ file or mapped BAM file |
| - | ref.fa | Input reference genome in FASTA format |
| - | working_dir | Specify working directory |

See [wiki](https://github.com/cytham/nanovar/wiki) for entire list of options.

### Output
| Output file | Comment |
| :--- | :--- |
| ${sample}.nanovar.pass.vcf | Final VCF filtered output file (1-based) |
| ${sample}.nanovar.pass.report.html | HTML report showing run summary and statistics |

For more information, see [wiki](https://github.com/cytham/nanovar/wiki).

### Full usage
```
usage: nanovar [options] [FASTQ/FASTA/BAM] [REFERENCE_GENOME] [WORK_DIRECTORY]

NanoVar is a neural network enhanced structural variant (SV) caller that handles low-depth long-read sequencing data.

positional arguments:
  [FASTQ/FASTA/BAM]     path to long reads or mapped BAM file.
                        Formats: fasta/fa/fa.gzip/fa.gz/fastq/fq/fq.gzip/fq.gz or .bam
  [reference_genome]    path to reference genome in FASTA. Genome indexes created
                        will overwrite indexes created by other aligners such as bwa.
  [work_directory]      path to work directory. Directory will be created
                        if it does not exist.

options:
  -h, --help            show this help message and exit
  --cnv hg38            also detects large genomic copy-number variations
                        using CytoCAD (e.g. loss/gain of whole chromosomes).
                        Only works with hg38 genome assembly. Please state 'hg38' [None]
  -x str, --data_type str
                        type of long-read data [ont]
                        ont - Oxford Nanopore Technologies
                        pacbio-clr - Pacific Biosciences CLR
                        pacbio-ccs - Pacific Biosciences CCS
  -f file, --filter_bed file
                        BED file with genomic regions to be excluded [None]
                        (e.g. telomeres and centromeres) Either specify name of in-built
                        reference genome filter (i.e. hg38, hg19, mm10) or provide full
                        path to own BED file.
  --annotate_ins str    enable annotation of INS with NanoINSight,
                        please specify species of sample [None]
                        Currently supported species are:
                        'human', 'mouse', and 'rattus'.
  -c int, --mincov int  minimum number of reads required to call a breakend [4]
  -l int, --minlen int  minimum length of SV to be detected [25]
  -p float, --splitpct float
                        minimum percentage of unmapped bases within a long read
                        to be considered as a split-read. 0.05<=p<=0.50 [0.05]
  -a int, --minalign int
                        minimum alignment length for single alignment reads [200]
  -b int, --buffer int  nucleotide length buffer for SV breakend clustering [50]
  -s float, --score float
                        score threshold for defining PASS/FAIL SVs in VCF [1.0]
                        Default score 1.0 was estimated from simulated analysis.
  --homo float          lower limit of a breakend read ratio to classify a homozygous state [0.75]
                        (i.e. Any breakend with homo<=ratio<=1.00 is classified as homozygous)
  --hetero float        lower limit of a breakend read ratio to classify a heterozygous state [0.35]
                        (i.e. Any breakend with hetero<=ratio<homo is classified as heterozygous)
  --debug               run in debug mode
  -v, --version         show version and exit
  -q, --quiet           hide verbose
  -t int, --threads int
                        number of available threads for use [1]
  --model path          specify path to custom-built model
  --mm path             specify path to 'minimap2' executable
  --st path             specify path to 'samtools' executable
  --ma path             specify path to 'mafft' executable for NanoINSight
  --rm path             specify path to 'RepeatMasker' executable for NanoINSight
```

### Operating system
* Linux (x86_64 architecture, tested in Ubuntu 14.04, 16.04, 18.04)  

### Installation
There are three ways to install NanoVar:
#### Option 1: Conda environment (Recommended)
```
conda create -n myenv -c bioconda python=3.11 samtools bedtools minimap2
conda activate myenv
pip install nanovar

or

conda create -n myenv -c bioconda python=3.11 nanovar
conda activate myenv
```
#### Option 2: PyPI (See dependencies below)
```
# Installing from PyPI requires own installation of dependencies, see below
pip install nanovar
```
#### Option 3: GitHub (See dependencies below)
```
# Installing from GitHub requires own installation of dependencies, see below
git clone https://github.com/cytham/nanovar.git 
cd nanovar 
pip install .
```
#### Installation of dependencies
* bedtools >=2.26.0
* samtools >=1.3.0
* minimap2 >=2.17
<!--
* makeblastdb and windowmasker
* hs-blastn ==0.0.5
-->

Please make sure each executable binary is in PATH.
##### 1. _bedtools_
Please visit [here](https://bedtools.readthedocs.io/en/latest/content/installation.html) for instructions to install.

##### 2. _samtools_
Please visit [here](http://www.htslib.org/download/) for instructions to install.

##### 3. _minimap2_
Please visit [here](https://github.com/lh3/minimap2) for instructions to install.
<!--
##### 4. _makeblastdb_ and _windowmasker_
```
# Download NCBI-BLAST v2.3.0+ from NCBI FTP server
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz

# Extract tar.gz
tar zxf ncbi-blast-2.3.0+-x64-linux.tar.gz

# Copy makeblastdb and windowmasker binaries to PATH (e.g. ~/bin)
cp ncbi-blast-2.3.0+/bin/makeblastdb ~/bin && cp ncbi-blast-2.3.0+/bin/windowmasker ~/bin
```
##### 5. _hs-blastn_
```
# Download and compile the 0.0.5 version
git clone https://github.com/chenying2016/queries.git
cd queries/hs-blastn-src/v0.0.5
make

# Copy hs-blastn binary to path (e.g. ~/bin)
cp hs-blastn ~/bin
```
If you encounter "isnan" error during compilation, please refer to [this](https://github.com/cytham/nanovar/issues/7#issuecomment-644546378).
-->

## Annotating INS variants with NanoINSight
NanoVar allows the concurrent repeat element annotation of INS variants using [NanoINSight](https://github.com/AsmaaSamyMohamedMahmoud/NanoINSight).

To run NanoINSight, simply add "--annotate_ins [species]" when running NanoVar.
```
nanovar -t 24 -f hg38 --annotate_ins human sample.bam ref.fa working_dir
```
To understand NanoINSight output files, please visit its repository [here](https://github.com/AsmaaSamyMohamedMahmoud/NanoINSight).

### Installation of NanoINSight dependencies

NanoINSight requires the installation of MAFFT and RepeatMasker. Please refer to [here](https://github.com/AsmaaSamyMohamedMahmoud/NanoINSight) for instructions on how to install them, or install them through Conda as shown below:

```
pip install nanoinsight
conda install -c bioconda mafft repeatmasker -y
```

Note: If encountered "numpy.dtype size changed" tensorflow error while running NanoVar, ensure numpy version is <2.0.0 (i.e. pip install numpy 1.26.4).
 
## Documentation
See [wiki](https://github.com/cytham/nanovar/wiki) for more information.

## Versioning
See [CHANGELOG](./CHANGELOG.txt)

## Citation
If you use NanoVar, please cite:

Tham, CY., Tirado-Magallanes, R., Goh, Y. et al. NanoVar: accurate characterization of patients’ genomic structural variants using low-depth nanopore sequencing. Genome Biol. 21, 56 (2020). https://doi.org/10.1186/s13059-020-01968-7


## Authors

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)
* **Roberto Tirado Magallanes** - [rtmag](https://github.com/rtmag)
* **Asmaa Samy** - [asmmahmoud](https://github.com/asmmahmoud)
* **Touati Benoukraf** - [benoukraflab](https://github.com/benoukraflab)

## License

This project is licensed under GNU General Public License - see [LICENSE.txt](./LICENSE.txt) for details.

## Simulation datasets and scripts used in the manuscript
SV simulation datasets used in the manuscript can be downloaded [here](https://doi.org/10.5281/zenodo.3569479). Scripts used for simulation dataset generation and tool performance comparison are available [here](./scripts).

Although NanoVar is provided with a universal model and threshold score, instructions required for building a custom neural-network model is available [here](https://github.com/cytham/nanovar/wiki/Model-training).

## Limitations
* The inaccurate basecalling of large homopolymer or low complexity DNA regions may result in the false determination of deletion SVs. We advise the use of up-to-date ONT basecallers such as [Dorado](https://github.com/nanoporetech/dorado) to minimize this possibility.

* For BND SVs, NanoVar is unable to calculate the actual number of SV-opposing reads (normal reads) at the novel adjacency as
 there are two breakends from distant locations. It is not clear whether the novel adjacency is derived from both or either
  breakends in cases of balanced and unbalanced variants, and therefore it is not possible to know which breakend location(s) to
   consider for counting normal reads. Currently, NanoVar approximates the normal read count by the minimum count from either 
   breakend location. Although this helps in capturing unbalanced BNDs, it might lead to some false positives.

## Tutorial for neural-network model training of NanoVar

### Introduction

The neural network model within the NanoVar SV calling pipeline functions to enrich true SV events based on read mapping characteristics and read support information. It provides a score for each SV which advise the confident level of its validity. This model was trained using simulated long-reads generated from a SV simulated genome. The model is applicable to real data assuming that the characteristics of real long-reads (e.g. error rate, indels) or genome where it came from (e.g. repetitive sequence proportions) are not vastly different from that of simulations. As the simulated long read charateristics were modeled after Oxford Nanopore (ONT) MinION long reads sequence from a human genome, we believe that the model would be applicable at least to ONT reads sequenced from a human genome. In situations where the model would be less effective, such as using reads from a novel long-read sequencer, we have provided this tutorial for users to train their own model based on the characteristics of their own reads and genome of their organism of interest. We have also provided methods to assess the model upon training completion and instructions on how to integrate the new model into NanoVar.


### Neural-network model training

#### Requirements

##### Softwares

* NanoVar (>= v1.2.8)
* RSVSim
* NanoSim
* minimap2

##### Python packages (Python 3)

* numpy
* matplotlib
* tensorflow (>= 2.0.0)

##### Data

* 500k to 2 million long-reads in FASTA/FASTQ (Could be less or more, depending on the user and organism)
* reference genome in FASTA (Please use only main contigs, remove any scaffold or unknown contigs)


#### Method overview

The first step is to create three genome FASTA files with the same simulated SVs but at different proportions. This allows the formation of homozygous/heterozygous SVs when long-reads are simulated from these genomes at specific proportions. For example, we simulate 10,000 SVs in Genome A, a subset of 8,000 SVs in Genome B, and a subset of 4,000 SVs in Genome C. The 4,000 SVs in Genome C are a subset of the 8,000 SVs in Genome B. Next, we simulate 2.5 million reads from Genome A, 2.5 million reads from Genome B and 5 million reads from Genome C. The combination of these 10 million reads creates 4,000 homozygous SV, 4,000 heterozygous SV, and 2,000 low-coverage SVs. A homozygous SV only has SV-supporting reads at their breakends while a heterozygous SV has both SV-supporting and SV-opposing reads at similar proportions. A low-coverage SV will have mostly SV-opposing reads with usually only one SV-supporting read. This combined read dataset will be used as the training set for the neural network model.

Method outline:
1. Simulate SV locations using a mock RSVSim run
2. Split CSV files to different SV proportions
3. Simulate SV genomes at different SV proportions using new CSV files
4. Simulate long-reads from the SV genomes to create the training dataset
5. Repeat steps 1-4 with a different seed to create the test dataset
6. Run NanoVar on the training and test dataset in --debug mode
7. Neural network training and assessment using tensorflow and keras
8. Integrate model into NanoVar for use in real data


#### 1. Simulate SV locations using a mock RSVSim run in R (https://bioconductor.org/packages/release/bioc/html/RSVSim.html)
```
# Load RSVSim library
library(RSVSim)

# Define your simulation settings
delno = 20000  # Number of deletion SVs
insno = 20000  # Combined number of viral insertion and transposition SVs
dupno = 5000  # Number of tandem duplication SVs
invno = 10000  # Number of inversion SVs

# Estimate SV sizes from real data
delSizes = estimateSVSizes(n=delno, minSize=50, maxSize=10000, default="deletions", hist=FALSE)
insSizes = estimateSVSizes(n=insno, minSize=50, maxSize=10000, default="insertions", hist=FALSE)
tandSizes = estimateSVSizes(n=dupno, minSize=50, maxSize=10000, default="tandemDuplications", hist=FALSE)
invSizes = estimateSVSizes(n=invno, minSize=50, maxSize=10000, default="inversions", hist=FALSE)

# Load own genome of reference
refgenome = readDNAStringSet("/path/to/own/genome.fa", format="fasta")

# Define output directory
out_dir = "/your/output/directory"

# Run simulation
sim = simulateSV(output=out_dir, genome=refgenome, del=delno, ins=insno, invs=invno, dups=dupno, sizeDels=delSizes, sizeIns=insSizes, sizeInvs=invSizes, sizeDups=tandSizes, maxDups=1, percCopiedIns=0, bpFlankSize=20, percSNPs=0.25, indelProb=0.5, maxIndelSize=5, bpSeqSize=200, seed=100, verbose=TRUE)

# Exit R, cd to directory containing the .csv files and proceed to step 2
```

#### 2. Split CSV files to different SV proportions
This step prepares the region files required for the actual SV genome simulation.
```
# Separate the insertions.csv to be used for transposition and viral insertion regions
tail -n +2 insertions.csv | head -n 10000 > viral_ins.csv
tail -n +10002 insertions.csv | cut -f 2-8 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$7"\t"$4"\t"$5"\t"$6}' > tpo.csv

# Use Python script "custom_ins_list_RSVSim.py" and "viral_genome.fa" FASTA to generate viral regions for insertion into
 regions specified in viral_ins.csv
python3 custom_ins_list_RSVSim.py total_viral_genome.fa viral_ins.csv 1 10000 1 ./

# Combine viral and tpo CSV files
cat viral_ins_regions1.csv tpo.csv > total_ins.csv

# Remove header from CSV files for deletion, inversion and duplication
tail -n +2 deletions.csv > del.csv 
tail -n +2 tandemDuplications.csv > dup.csv 
tail -n +2 inversions.csv > inv.csv 

```

#### 3. Simulate SV genomes at different SV proportions using new CSV files


#### 4. Simulate long-reads from the SV genomes to create the training dataset


#### 5. Repeat steps 1-4 with a different seed to create the test dataset


#### 6. Run NanoVar on the training and test dataset in --debug mode


#### 7. Neural network training and assessment using tensorflow and keras


#### 8. Integrate model into NanoVar for use in real data

"""
Functions for parsing and verifying input parameters.

Copyright (C) 2021 Tham Cheng Yong

This file is part of NanoVar.

NanoVar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

NanoVar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with NanoVar.  If not, see <https://www.gnu.org/licenses/>.
"""


import sys
import argparse
from nanovar import __version__


# Parse input
def input_parser(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="NanoVar is a neural network enhanced structural variant (SV) caller that \
handles low-depth long-read sequencing data.",
                                     formatter_class=argparse.RawTextHelpFormatter, usage=msg())  # RawDescriptionHelpFormatter)

    def restrict_float(f):
        f = float(f)
        if f < 0.05 or f > 0.50:
            raise argparse.ArgumentTypeError("%r not in range [0.05, 0.50]" % (f,))
        return f

    parser.add_argument("input", type=str,
                        metavar="[FASTQ/FASTA/BAM]",
                        help="""path to long reads or mapped BAM file.
Formats: fasta/fa/fa.gzip/fa.gz/fastq/fq/fq.gzip/fq.gz or .bam""")

    parser.add_argument("ref", type=str,
                        metavar="[reference_genome]",
                        help="""path to reference genome in FASTA. Genome indexes created 
will overwrite indexes created by other aligners such as bwa.""")

    parser.add_argument("dir", type=str,
                        metavar="[work_directory]",
                        help="""path to work directory. Directory will be created 
if it does not exist.""")

    parser.add_argument("--cnv", type=str, metavar="hg38",
                        default=None,
                        help="""also detects large genomic copy-number variations 
using CytoCAD (e.g. loss/gain of whole chromosomes). 
Only works with hg38 genome assembly. Please state 'hg38' [None]""")

    parser.add_argument("-x", "--data_type", type=str, metavar="str",
                        default='ont',
                        help="""type of long-read data [ont]
ont - Oxford Nanopore Technologies
pacbio-clr - Pacific Biosciences CLR
pacbio-ccs - Pacific Biosciences CCS""")

    parser.add_argument("-f", "--filter_bed", type=str, metavar="file",
                        help="""BED file with genomic regions to be excluded [None]
(e.g. telomeres and centromeres) Either specify name of in-built 
reference genome filter (i.e. hg38, hg19, mm10) or provide full 
path to own BED file.""")

    parser.add_argument("--annotate_ins", type=str, metavar="str",
                        default=None,
                        help="""enable annotation of INS with NanoINSight, 
please specify species of sample [None]
Currently supported species are:
'human', 'mouse', and 'rattus'.
""")

    parser.add_argument("-c", "--mincov", type=int, metavar="int",
                        default=2,
                        help="minimum number of reads required to call a breakend [2]")

    parser.add_argument("-l", "--minlen", type=int, metavar="int",
                        default=25,
                        help="minimum length of SV to be detected [25]")

    parser.add_argument("-p", "--splitpct", type=restrict_float, metavar="float",
                        default=0.05,
                        help="""minimum percentage of unmapped bases within a long read 
to be considered as a split-read. 0.05<=p<=0.50 [0.05]""")

    parser.add_argument("-a", "--minalign", type=int, metavar="int",
                        default=200,
                        help="minimum alignment length for single alignment reads [200]")

    parser.add_argument("-b", "--buffer", type=int, metavar="int",
                        default=50,
                        help="nucleotide length buffer for SV breakend clustering [50]")

    parser.add_argument("-s", "--score", type=float, metavar="float",
                        default=1.0,
                        help="""score threshold for defining PASS/FAIL SVs in VCF [1.0]
Default score 1.0 was estimated from simulated analysis. """)

    parser.add_argument("--homo", type=float, metavar="float",
                        default=0.75,
                        help="""lower limit of a breakend read ratio to classify a homozygous state [0.75]
(i.e. Any breakend with homo<=ratio<=1.00 is classified as homozygous)""")

    parser.add_argument("--hetero", type=float, metavar="float",
                        default=0.35,
                        help="""lower limit of a breakend read ratio to classify a heterozygous state [0.35]
(i.e. Any breakend with hetero<=ratio<homo is classified as heterozygous)""")

    parser.add_argument("--debug", action='store_true',
                        help="run in debug mode")

    # parser.add_argument("--force", action='store_true',
    #                     help="run full pipeline (i.e. do not skip index generation)")

    parser.add_argument("-v", "--version", action='version',
                        version=__version__,
                        help="show version and exit")

    parser.add_argument("-q", "--quiet", action='store_true',
                        help="hide verbose")

    parser.add_argument("-t", "--threads", type=int, metavar="int",
                        default=1,
                        help="number of available threads for use [1]")

    parser.add_argument("--model", type=str, metavar="path",
                        help="specify path to custom-built model")

    parser.add_argument("--mm", type=str, metavar="path",
                        help="specify path to 'minimap2' executable")

    parser.add_argument("--st", type=str, metavar="path",
                        help="specify path to 'samtools' executable")

    parser.add_argument("--ma", type=str, metavar="path",
                        help="specify path to 'mafft' executable for NanoINSight")

    parser.add_argument("--rm", type=str, metavar="path",
                        help="specify path to 'RepeatMasker' executable for NanoINSight")

    # parser.add_argument("--mdb", type=str, metavar="path",
    #                     help="specify path to 'makeblastdb' executable")

    # parser.add_argument("--wmk", type=str, metavar="path",
    #                     help="specify path to 'windowmasker' executable")

    # parser.add_argument("--hsb", type=str, metavar="path",
    #                     help="specify path to 'hs-blastn' executable")

    parser.add_argument("--pickle", action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument("--archivefasta", action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument("--blastout", type=str, metavar="path",
                        help=argparse.SUPPRESS)

    args = parser.parse_args(args)
    return args


# Check if fasta/fastq is gzip compressed
def gzip_check(path):
    twobytes = b'\x1f\x8b'
    with open(path, 'rb') as f:
        return f.read(2) == twobytes


# Custom usage message
def msg():
    return "nanovar [options] [FASTQ/FASTA/BAM] [REFERENCE_GENOME] [WORK_DIRECTORY]"

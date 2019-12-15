"""
Functions for parsing and verifying input parameters.

Copyright (C) 2019 Tham Cheng Yong

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
    parser = argparse.ArgumentParser(description="""NanoVar is a neural-network-based structural variant (SV) caller that utilizes 
low-depth long-read sequencing data.""",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    def restrict_float(f):
        f = float(f)
        if f < 0.05 or f > 0.50:
            raise argparse.ArgumentTypeError("%r not in range [0.05, 0.50]" % (f,))
        return f

    parser.add_argument("reads", type=str,
                        metavar="[long_reads]",
                        help="Path to long reads. \
                             Formats: fasta/fa/fa.gzip/fa.gz \
                             fastq/fq/fq.gzip/fq.gz")

    parser.add_argument("ref", type=str,
                        metavar="[reference_genome]",
                        help="Path to reference genome in FASTA. \
                             Genome indexes created will overwrite \
                             indexes created by other aligners (e.g. bwa)")

    parser.add_argument("dir", type=str,
                        metavar="[working_directory]",
                        help="Path to working directory. \
                             Directory will be created if it does not exist")

    parser.add_argument("-f", "--filter_bed", type=str,
                        help="BED file with genomic regions to be excluded. \
                             (e.g. telomeres and centromeres). \
                             Either specify name of in-built reference genome filter \
                             (i.e. hg38, hg19, mm10) or provide FULL path to \
                             own BED file. [None]")

    parser.add_argument("-c", "--mincov", type=int,
                        default=1,
                        help="minimum number of reads required to call a breakend [1]")

    parser.add_argument("-l", "--minlen", type=int,
                        default=25,
                        help="minimum length of SV to be detected. [25]")

    parser.add_argument("-p", "--splitpct", type=restrict_float,
                        default=0.05,
                        help="minimum percentage of unmapped bases within a long read \
                             to be considered as a split-read. 0.05<=p<=0.50 [0.05]")

    parser.add_argument("-a", "--minalign", type=int,
                        default=200,
                        help="minimum alignment length for single alignment reads. [200]")

    parser.add_argument("-b", "--buffer", type=int,
                        default=50,
                        help="nucleotide length buffer for SV breakend clustering. [50]")

    parser.add_argument("-s", "--score", type=float,
                        default=1.0,
                        help="score threshold for defining PASS/FAIL SVs in VCF. \
                             Default score 1.0 was estimated from simulated analysis. [1.0]")

    parser.add_argument("--homo", type=float,
                        default=0.75,
                        help="Lower limit of a breakend read ratio to classify a homozygous state [0.75]")

    parser.add_argument("--hetero", type=float,
                        default=0.35,
                        help="Lower limit of a breakend read ratio to classify a heterozygous state [0.35]")

    parser.add_argument("-v", "--version", action='version',
                        version=__version__,
                        help="show version and exit")

    parser.add_argument("-q", "--quiet", action='store_true',
                        help="hide verbose")

    parser.add_argument("-t", "--threads", type=int, choices=range(1, 55), metavar="[1-54]",
                        default=2,
                        help="number of available threads for use, max=54 [2]")

    parser.add_argument("--force", action='store_true',
                        help="run full pipeline, do not skip any redundant steps (e.g. index generation)")

    parser.add_argument("--mdb", type=str,
                        help="Specify path to 'makeblastdb' executable")

    parser.add_argument("--wmk", type=str,
                        help="Specify path to 'windowmasker' executable")

    parser.add_argument("--hsb", type=str,
                        help="Specify path to 'hs-blastn' executable")

    args = parser.parse_args(args)
    return args


# Check if fasta/fastq is gzip compressed
def gzip_check(path):
    twobytes = b'\x1f\x8b'
    with open(path, 'rb') as f:
        return f.read(2) == twobytes

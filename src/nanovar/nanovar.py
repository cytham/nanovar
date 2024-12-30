#!/usr/bin/env python3

"""
NanoVar

This is the main executable file of the program NanoVar.

Copyright (C) 2021 Tham Cheng Yong

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


__author__ = 'CY Tham'
import os
import sys
import time
import pysam
import shutil
import logging
import nanovar
import random
import pickle
import threading
from datetime import datetime
from nanovar import __version__, input_parser, gzip_check, fastx_valid, bed_valid, check_exe, align_mm
#from progress.spinner import Spinner
#import nanoinsight


def main():
    # Parse arguments
    args = input_parser()
    file_path = args.input
    ref_path = args.ref
    ref_name = os.path.basename(ref_path).rsplit('.', 1)[0]
    wk_dir = args.dir
    data_type = args.data_type
    genome_filter = args.filter_bed
    cnv = args.cnv
    annotate_ins = args.annotate_ins
    minlen = args.minlen
    splitpct = args.splitpct
    minalign = args.minalign
    mincov = args.mincov
    buff = args.buffer
    score_threshold = args.score
    homo_t = args.homo
    het_t = args.hetero
    threads = args.threads
    debug = args.debug
    quiet = args.quiet
    # force = args.force
    model_path = args.model
    filter_bed_dir = os.path.join(os.path.dirname(nanovar.__file__), 'gaps')
    ref_dir = os.path.join(os.path.dirname(nanovar.__file__), 'ref')
    mm = args.mm
    st = args.st
    ma = args.ma
    rm = args.rm
    # mdb = args.mdb
    # wmk = args.wmk
    # hsb = args.hsb
    pickle_bool = args.pickle
    archivefasta = args.archivefasta
    # blastout = args.blastout
    blastout = 'blastoff'
    nv_cmd = ' '.join(sys.argv)

    # Check model file
    if model_path is None:
        if data_type == 'ont':
            model_path = os.path.join(os.path.dirname(nanovar.__file__), 'model',
                                      'ANN.E100B400L3N12-5D0.4-0.2SGDsee11_het_gup_v1.h5')
        elif data_type == 'pacbio-clr':
            model_path = os.path.join(os.path.dirname(nanovar.__file__), 'model',
                                      'ANN.E100B400L3N12-5D0.4-0.2SGDsee11_het_clr_v1.h5')
        elif data_type == 'pacbio-ccs':
            model_path = os.path.join(os.path.dirname(nanovar.__file__), 'model',
                                      'ANN.E100B400L3N12-5D0.4-0.2SGDsee11_het_ccs_v1.h5')
        else:
            logging.critical("Error: Invalid data type given '%s'" % data_type)
            raise Exception("Error: Invalid data type given '%s'" % data_type)
    else:
        if data_type not in ['ont', 'pacbio-clr', 'pacbio-ccs']:
            logging.info("Invalid data type given '%s', but irrelevant due to custom-built model use." % data_type)

    # Check homo_t > het_t
    if homo_t <= het_t:
        raise Exception("Error: --homo threshold %s is less than or equal to --hetero threshold %s" % (str(homo_t), str(het_t)))

    # Check for required executables
    mm = check_exe(mm, 'minimap2')
    st = check_exe(st, 'samtools')
    # mdb = check_exe(mdb, 'makeblastdb')
    # wmk = check_exe(wmk, 'windowmasker')
    # hsb = check_exe(hsb, 'hs-blastn')

    # Pre-check if annotating INS with NanoINSight
    if annotate_ins:
        try:
            import nanoinsight
        except ImportError:
            raise Exception("Error: nanoinsight module not found. Please install by 'pip install nanoinsight'.")
        species = annotate_ins
        nanoinsight.check_args(species)
        mafft_exe = nanoinsight.check_exe(ma, 'mafft')
        repmask_exe = nanoinsight.check_exe(rm, 'RepeatMasker')

    # If additional CNV mode
    hg38_size_dict = {}
    if cnv:
        raise Exception("Error: CNV module presently not available on this version. Please remove '--cnv' option. Sorry, we are actively improving.")
        # # Check CytoCAD
        # try:
        #     import cytocad
        #     data_dir = os.path.join(os.path.dirname(cytocad.__file__), 'data')
        #     hg38_size_path = os.path.join(data_dir, 'hg38_sizes_main.bed')
        #     with open(hg38_size_path) as f:
        #         for bed in f:
        #             line = bed.split('\t')
        #             hg38_size_dict[line[0]] = int(line[2])
        # except ImportError:
        #     raise Exception("Error: cytocad module not found, please install CytoCAD - https://github.com/cytham/cytocad")
        # # Check tagore and rsvg
        # _ = check_exe(None, 'tagore')
        # _ = check_exe(None, 'rsvg-convert')
        # # Check hg38 build input
        # if cnv != 'hg38':
        #     raise Exception("Error: Reference genome build %s is not recognised. CytoCAD only supports build hg38." % cnv)

    # Observe verbosity
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    # Assign threads
    threads_mm = threads

    # Print initiation message
    now = datetime.now()
    now_str = now.strftime("[%d/%m/%Y %H:%M:%S]")
    print(now_str, "- NanoVar-v%s started" % __version__)
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Checking integrity of input files - ', end='', flush=True)

    # Setup working directory
    if not os.path.exists(wk_dir):
        os.makedirs(wk_dir)
    if not os.path.exists(os.path.join(wk_dir, 'fig')):
        os.makedirs(os.path.join(wk_dir, 'fig'))

    # Setup up logging
    log_file = os.path.join(wk_dir, 'NanoVar-{:%d%m%y-%H%M}.log'.format(datetime.now()))
    logging.basicConfig(filename=log_file, level=logging.DEBUG, format='[%(asctime)s] - %(levelname)s - %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')
    logging.info('Initialize NanoVar log file')
    logging.info('Version: NanoVar-%s' % __version__)
    logging.info('Command: %s' % nv_cmd)

    # Detect read or map file
    filename = os.path.basename(file_path)
    read_suffix = ['.fa', '.fq', '.fasta', '.fastq', '.fa.gzip', '.fq.gzip', '.fa.gz', '.fq.gz', '.fasta.gz', '.fastq.gz']
    bam_suffix = '.bam'
    cram_suffix = '.cram'
    contig_list = []
    if any(filename.lower().endswith(s) for s in read_suffix):
        input_name = os.path.basename(file_path).rsplit('.f', 1)[0]
        input_type = 'raw'
        # Test gzip compression and validates read file
        if gzip_check(file_path):
            # read_para = "<(zcat " + file_path + ")"
            fastx_check = fastx_valid(file_path, "gz")
        else:
            # read_para = "<(cat " + file_path + ")"
            fastx_check = fastx_valid(file_path, "txt")
        if fastx_check[0] == "Fail":
            logging.critical("Error: Input FASTQ/FASTA file is corrupted around line %s +/- 4" % str(fastx_check[1]))
            raise Exception("Error: Input FASTQ/FASTA file is corrupted around line %s +/- 4" % str(fastx_check[1]))
        elif fastx_check[0] == "Pass":
            logging.debug("Input FASTQ/FASTA file passed")
    elif filename.lower().endswith(bam_suffix):
        save = pysam.set_verbosity(0)  # Suppress BAM index missing warning
        sam = pysam.AlignmentFile(file_path, "rb")
        pysam.set_verbosity(save)  # Revert verbosity level
        try:
            assert sam.is_bam, "Error: Input BAM file is not a BAM file."
            input_name = os.path.basename(file_path).rsplit('.bam', 1)[0]
            input_type = 'bam'
            fastx_check = []
            # Get BAM contigs from header
            header = sam.header.to_dict()
            for h in header['SQ']:
                contig_list.append(h['SN'])
        except AssertionError:
            logging.critical("Error: Input BAM file is not a BAM file.")
            raise Exception("Error: Input BAM file is not a BAM file.")
    elif filename.lower().endswith(cram_suffix):
        save = pysam.set_verbosity(0)  # Suppress index missing warning
        sam = pysam.AlignmentFile(file_path, "rc", reference_filename=ref_path)
        pysam.set_verbosity(save)  # Revert verbosity level
        try:
            assert sam.is_cram, "Error: Input CRAM file is not a CRAM file."
            input_name = os.path.basename(file_path).rsplit('.cram', 1)[0]
            input_type = 'cram'
            fastx_check = []
            # Get CRAM contigs from header
            header = sam.header.to_dict()
            for h in header['SQ']:
                contig_list.append(h['SN'])
        except AssertionError:
            logging.critical("Error: Input CRAM file is not a CRAM file.")
            raise Exception("Error: Input CRAM file is not a CRAM file.")
    else:
        logging.critical("Error: Input file is not recognised, please ensure file suffix has '.fa' or '.fq' or '.bam' or '.cram'")
        raise Exception("Error: Input file is not recognised, please ensure file suffix has '.fa' or '.fq' or '.bam' or '.cram'")

    if gzip_check(ref_path):
        logging.critical("Error: Input reference file is gzipped, please unzip it")
        raise Exception("Error: Input reference file is gzipped, please unzip it")

    # Logging config info
    logging.info('Input file: %s' % file_path)
    logging.info('Read type: %s' % data_type)
    logging.info('Reference genome: %s' % ref_path)
    logging.info('Working directory: %s' % wk_dir)
    logging.info('minimap2 path: %s' % shutil.which(mm))
    logging.info('samtools path: %s' % shutil.which(st))
    logging.info('Model: %s' % model_path)
    logging.info('Filter file: %s' % genome_filter)
    logging.info('Minimum number of reads for calling a breakend: %s' % str(mincov))
    logging.info('Minimum SV len: %s' % str(minlen))
    logging.info('Mapping percent for split-read: %s' % str(splitpct))
    logging.info('Length buffer for clustering: %s' % str(buff))
    logging.info('Score threshold: %s' % str(score_threshold))
    logging.info('Homozygous read ratio threshold: %s' % str(homo_t))
    logging.info('Heterozygous read ratio threshold: %s' % str(het_t))
    logging.info('Number of threads: %s\n' % str(threads))
    if input_type == 'raw':
        logging.info('Total number of reads in FASTQ/FASTA: %s\n' % str(fastx_check[1]))
    elif input_type == 'bam' or input_type == 'cram':
        logging.info('Total number of reads in FASTQ/FASTA: -\n')
    logging.info('NanoVar started')

    from Bio import SeqIO
    from collections import OrderedDict
    # Process reference genome
    contig_len_dict = OrderedDict()
    total_gsize = 0
    for seq_record in SeqIO.parse(ref_path, "fasta"):
        contig_len_dict[seq_record.id] = len(seq_record)
        total_gsize += len(seq_record)

    # Check if contig sizes from BAM/CRAM is same as hg38 if CNV mode
    if cnv:
        for chrom in hg38_size_dict:
            if chrom in contig_len_dict:
                if contig_len_dict[chrom] != hg38_size_dict[chrom]:
                    raise Exception("Error: Contig %s in BAM/CRAM (%i) has different length in hg38 genome assembly (%i)" %
                                    (chrom, contig_len_dict[chrom], hg38_size_dict[chrom]))
            else:
                raise Exception("Error: hg38 contig %s is absent in BAM/CRAM" % chrom)

    # Check BAM/CRAM contigs in reference genome
    if input_type == 'bam' or input_type == 'cram':
        for c in contig_list:
            if c not in contig_len_dict:
                logging.critical("Error: Contig %s in BAM/CRAM is absent in reference genome" % c)
                raise Exception("Error: Contig %s in BAM/CRAM is absent in reference genome" % c)

    # Check contig id for invalid symbols
    contig_omit = checkcontignames(contig_len_dict)

    # Validate filter BED file
    if genome_filter is not None:
        if genome_filter == 'hg38':
            filter_path = os.path.join(filter_bed_dir, genome_filter + '_curated_filter_main.bed')
        elif genome_filter in ('hg19', 'mm10'):
            filter_path = os.path.join(filter_bed_dir, genome_filter + '_filter-B400.bed')
        else:
            filter_path = genome_filter
        if os.path.isfile(filter_path):
            if bed_valid(filter_path, contig_len_dict):
                logging.debug("Genome filter BED passed")
        else:
            logging.critical("Error: Genome filter BED %s is not found" % filter_path)
            raise Exception("Error: Genome filter BED %s is not found" % filter_path)
    else:
        filter_path = genome_filter

    print('Pass')

    # Update progress
    if len(contig_omit) > 0:
        print('Warning: The following contig id(s) contains invalid symbols [ ~ : - ] %s. Reads mapping to these contig(s) '
              'will be ignored.' % ', '.join(['>' + x for x in contig_omit]))
        logging.warning('Warning: The following contig id(s) contains invalid symbols [ ~ : - ] %s. Reads mapping to these '
                        'contig(s) will be ignored.' % ', '.join(['>' + x for x in contig_omit]))

    # Aligning using minimap2
    if input_type == 'raw':
        logging.info('Read alignment using minimap2')
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Read mapping using minimap2 - ', end='', flush=True)
        mma = align_mm(ref_path, file_path, wk_dir, input_name, ref_name, threads_mm, mm, data_type, st)
        bam_path = mma[1]
        save = pysam.set_verbosity(0)  # Suppress index missing warning
        sam = pysam.AlignmentFile(bam_path, "rb")
        pysam.set_verbosity(save)  # Revert verbosity level
        print('Done')
    elif input_type == 'bam' or input_type == 'cram':
        logging.info('Input BAM/CRAM file, skipping minimap2 alignment')
        mma = ['-', '']
        bam_path = file_path

    # Update progress
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Analyzing read alignments and detecting SVs - ', end='', flush=True)
    
    # Parse and detect SVs
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    from nanovar.nv_characterize import VariantDetect
    logging.info('Parsing BAM and detecting SVs')
    run = VariantDetect(wk_dir, bam_path, splitpct, minalign, filter_path, minlen, buff, model_path,
                        total_gsize, contig_len_dict, score_threshold, file_path, input_name, ref_path, ref_name, mma[0],
                        mincov, homo_t, het_t, debug, contig_omit, cnv, nv_cmd, sam)
    run.bam_parse_detect()
    run.coverage_stats()
    print('Done')

    logging.info('Total number of mapped reads: %s\n' % str(run.maps))
    
    # Update progress
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Clustering SV breakends - ', end='', flush=True)

    # SV breakend clustering and extracting INS and INV SVs
    run.cluster_nn2()
    print('Done')

    # Analyze TE
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Correcting DUP and detecting TE - ', end='', flush=True)
    pysam.faidx(os.path.join(wk_dir, 'temp1.fa'))
    run.dup_te_detect(ref_dir, threads, mm, st, data_type) #

    # Write parse2 and cluster intermediate files
    run.write2file(add_out=[])
    print('Done')

    # Pickle objects
    if pickle_bool:
        #saveobjects = (wk_dir, bam_path, splitpct, minalign, filter_path, minlen, buff, model_path,
        #               total_gsize, contig_len_dict, score_threshold, file_path, input_name, ref_path, ref_name, mma[0],
        #               mincov, homo_t, het_t, debug, contig_omit, cnv, run.rlendict, run.total_subdata, run.total_out,
        #               run.out_nn, sub_run.total_out)
        run.sam = None # Remove pysam instance 
        with open(os.path.join(wk_dir, "nv_run.pkl"), "wb") as f:
            pickle.dump(run, f)
    
    # Update progress
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Generating VCF files and report - ', end='', flush=True)
    
    # Generating VCF and HTML report
    run.vcf_report()
    print('Done')

    # Annotate insertions
    if annotate_ins:
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Annotating INS with NanoINSight - ', end='', flush=True)
        logging.info('Annotating INS with NanoINSight')
        vcf = os.path.join(wk_dir, '%s.nanovar.pass.vcf' % input_name)
        sv_sup = os.path.join(wk_dir, 'sv_support_reads.tsv')
        ins_seq = os.path.join(wk_dir, 'ins_seq.fa')
        id_seq, fasta_dir = nanoinsight.create_fa(vcf, wk_dir, sv_sup, ins_seq)
        con_fasta, threads_per_job = nanoinsight.create_cons(vcf, wk_dir, fasta_dir, id_seq, threads, mafft_exe, batch_size=100, num_parallel_workers=5)
        nanoinsight.rep_annote(wk_dir, con_fasta, threads_per_job, species, repmask_exe)
        print('Done')

    # Output SV-supporting BAM
    if args.sv_bam_out:
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"), '- Creating SV-supporting BAM - ', end='', flush=True)
        logging.info('Creating SV-supporting BAM')
        from .nv_supp_bam import create_sv_supp_bam
        vcf = os.path.join(wk_dir, '%s.nanovar.pass.vcf' % input_name)
        sv_sup = os.path.join(wk_dir, 'sv_support_reads.tsv')
        create_sv_supp_bam(vcf, sv_sup, bam_path, wk_dir, input_type)
        print('Done')
    
    # Delete temporary fasta file
    if not archivefasta:
        f = os.path.join(wk_dir, 'temp1.fa')
        os.remove(f)
        f = os.path.join(wk_dir, 'temp1.fa.fai')
        os.remove(f)
    
    if not debug:
        f = os.path.join(wk_dir, 'temp2.fa')
        os.remove(f)
        f = os.path.join(wk_dir, '%s.nanovar.total.vcf' % input_name)
        os.remove(f)
        f = os.path.join(wk_dir, 'genome.sizes')
        os.remove(f)
    
    logging.info('NanoVar ended')
    now = datetime.now()
    now_str = now.strftime("[%d/%m/%Y %H:%M:%S]")
    print(now_str, "- NanoVar ended")

# Check contig name
def checkcontignames(contig_len_dict):
    contig_omit = {}
    for contig in contig_len_dict:
        if "~" in contig or ":" in contig or "-" in contig:
            contig_omit[contig] = [0, int(contig_len_dict[contig])]
    return contig_omit

if __name__ == "__main__":
    main()

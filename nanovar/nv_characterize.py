"""
Method for SV characterization

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


import os
import logging
import random
from nanovar.nv_bam_parser import bam_parse
# from nanovar.nv_detect_algo import sv_detect
# from nanovar.nv_parser import entry_parser, align_info, breakpoint_parser
from nanovar.nv_cluster import sv_cluster
from nanovar.nv_nn import inference, svread_ovl
from nanovar.nv_cov_upper import ovl_upper
from nanovar.nv_vcf import create_vcf
from nanovar.nv_report import create_report
from nanovar.nv_dup_te_detect import dup_te_analyzer
from nanovar.nv_alt_seq import get_alt_seq
#from cytocad.change_detection import cad
#from cytocad.ideogram import tagore_wrapper


class VariantDetect:

    def __init__(self, wk_dir, bam, splitpct, minalign, filter_path, minlen, buff, model_path, total_gsize,
                 contig_len_dict, thres, read_path, read_name, ref_path, ref_name, map_cmd, mincov, homo_t, het_t, debug,
                 contig_omit, cnv, nv_cmd, sam):
        self.dir = wk_dir
        self.bam = bam
        self.sam = sam
        self.splitpct = splitpct
        self.minalign = minalign
        self.filter = filter_path
        self.contig_omit = contig_omit
        self.minlen = minlen
        self.buff = buff
        self.model = model_path
        self.gsize = total_gsize
        self.contig = contig_len_dict
        self.thres = thres
        self.rpath = read_path
        self.rname = read_name
        self.refpath = ref_path
        self.refname = ref_name
        self.mapcmd = map_cmd
        self.mincov = mincov
        self.homo_t = homo_t
        self.het_t = het_t
        self.debug = debug
        self.cnv = cnv
        self.basecov, self.maxovl, self.depth, self.maxovl3 = 0, 0, 0, 0
        self.total_out, self.total_subdata, self.out_nn, self.ins_out, self.out_rest, self.detect_out, self.beddata = [], [], [], [], [], [], []
        self.rlendict, self.parse_dict, self.index2te = {}, {}, {}
        # HTML SV table entry limit
        self.num_limit = 1000
        # HTML SV table SV ratio limit
        self.ratio_limit = 1
        self.maps = 0
        self.seed = 0
        self.seed2 = 1
        self.nv_cmd = nv_cmd

    def bam_parse_detect(self):
        random.seed(1)
        self.total_subdata, self.total_out, self.basecov, self.parse_dict, self.rlendict, self.maps, self.detect_out, self.seed, self.beddata \
            = bam_parse(self.sam, self.minlen, self.splitpct, self.minalign, self.dir, self.filter, self.contig_omit)
        writer(os.path.join(self.dir, 'subdata.tsv'), self.total_subdata, self.debug)
        writer(os.path.join(self.dir, 'detect.tsv'), self.detect_out, self.debug)
        writer(os.path.join(self.dir, 'parse1.tsv'), self.total_out, self.debug)

    def coverage_stats(self):
        # Obtaining upper cov limit and depth of coverage
        self.maxovl, self.depth, self.maxovl3 = ovl_upper(self.gsize, self.contig, self.basecov, self.beddata, self.dir)
        # CNV detection using cytocad
        if self.cnv:
            cnv_out, tag = cad(self.total_subdata, self.depth, self.maxovl3, self.rname, ref_build=self.cnv, cov_plots=True,
                               wk_dir=self.dir)
            tagore_wrapper(tag, self.rname, self.dir, self.cnv, 'png')
            cnv_out_path = os.path.join(self.dir, self.rname + ".CNV.bed")
            outwrite = open(cnv_out_path, 'w')
            _ = outwrite.write('\n'.join(cnv_out))
            outwrite.close()
            logging.info("CNV analysis completed.")
        # Report statistics of coverage
        logging.info("Genome size: %s bases" % str(self.gsize))
        logging.info("Mapped bases: %s bases" % str(self.basecov))
        logging.info("Depth of coverage: %sx" % str(self.depth))
        if self.depth < 4:
            logging.warning("Sequencing depth is less than 4x, output may not be comprehensive")
        logging.info("Read overlap upper limit: %s\n" % str(self.maxovl))

    def cluster_nn2(self):
        logging.info("Clustering SV breakends")
        # cluster_out, _ = sv_cluster(self.beddata, self.total_out, self.buff, self.maxovl, self.mincov,
        #                            self.contig, True, self.seed2)
        # Removed maxovl filter
        cluster_out, _ = sv_cluster(self.beddata, self.total_out, self.buff, 1000000, self.mincov,
                                    self.contig, True, self.seed2)
        logging.info("Neural network inference")
        if cluster_out:
            self.out_nn = inference(cluster_out, self.total_out, self.model)
        else:
            self.out_nn = []
        # Generate sv_overlap file
        #svread_ovl(self.dir, self.out_nn)

    # def cluster_extract(self):
    #     logging.info("Clustering SV breakends")
    #     cluster_parse = sv_cluster(self.total_subdata, self.total_out, self.buff, self.maxovl, self.mincov,
    #                                self.contig, False, self.seed2)
    #     logging.info("Filtering INS and INV SVs")
    #     total_qnames = []
    #     for line in cluster_parse:
    #         svtype = line.split('~')[1]
    #         qnames = []
    #         if svtype in ['bp_Nov_Ins', 'Nov_Ins']:  # Some INS are actually DUP
    #             for r in cluster_parse[line]:
    #                 qnames.append(r.rsplit('~', 1)[0])
    #                 self.ins_out.append(self.parse_dict[r])
    #             total_qnames.extend(qnames)
    #         elif svtype == 'Del':
    #             for r in cluster_parse[line]:
    #                 self.ins_out.append(self.parse_dict[r])
    #         else:
    #             for r in cluster_parse[line]:
    #                 qnames.append(r.rsplit('~', 1)[0])
    #             total_qnames.extend(qnames)
    #     qnames_dict = {x: 1 for x in set(total_qnames)}
    #     self.fasta_extract(qnames_dict)

    # def cluster_nn(self, add_out):
    #     logging.info("Re-clustering INS/INV SVs and merging")
    #     # Merge old ins with new from hsblastn
    #     sub_out = self.ins_out + add_out
    #     cluster_out_ins, _ = sv_cluster(self.total_subdata, sub_out, self.buff, self.maxovl, self.mincov, self.contig, True,
    #                                     self.seed2)
    #     new_cluster_out = self.out_rest + cluster_out_ins
    #     logging.info("Neural network inference")
    #     new_total_out = self.total_out + add_out
    #     if new_cluster_out:
    #         self.out_nn = inference(new_cluster_out, new_total_out, self.model)
    #     else:
    #         self.out_nn = []
    #     # Generate sv_overlap file
    #     svread_ovl(self.dir, self.out_nn)

    # def parse_detect_hsb(self):
    #     # Make gap dictionary
    #     gapdict = makegapdict(self.filter, self.contig_omit)
    #     temp1 = []
    #     chromocollect = []
    #     ovlt = 0.9
    #     sig_index = [0, 2, 4, 6, 8]
    #     data = open(self.bam, 'r').read().splitlines()
    #     data.append('null\tnull\tnull\tnull\tnull\tnull')
    #     nlines = len(data) - 1
    #     for i in range(nlines):
    #         if data[i].split('\t')[0] == data[i + 1].split('\t')[0]:  # Grouping alignments by read name
    #             temp1.append(align_info(data[i], self.rlendict))
    #             if data[i].split('\t')[1].strip() not in chromocollect:
    #                 chromocollect.append(data[i].split('\t')[1].strip())
    #         else:
    #             temp1.append(align_info(data[i], self.rlendict))
    #             if data[i].split('\t')[1].strip() not in chromocollect:
    #                 chromocollect.append(data[i].split('\t')[1].strip())
    #             self.seed += 1
    #             # Parse entries and correct overlap alignments
    #             subdata = entry_parser(temp1, chromocollect, ovlt)
    #             # SV detection
    #             out1, out2 = sv_detect(subdata, self.splitpct, self.minalign, gapdict)
    #             if out2 == '':
    #                 pass
    #             else:
    #                 # Parse breakpoints
    #                 final = breakpoint_parser(out2, self.minlen, sig_index, self.seed, 'hsb')
    #                 self.total_out.extend(final)
    #             temp1 = []
    #             chromocollect = []
    #     if not self.debug:  # Remove blast table if not debug mode
    #         os.remove(self.bam)
    
    def dup_te_detect(self, ref_dir, threads, mm, st, data_type):
        logging.info("Detecting DUP and TE")
        self.index2te, self.out_nn = dup_te_analyzer(self.dir, self.out_nn, self.total_out, self.thres, ref_dir, self.refpath, mm, threads, data_type, st, self.debug)
    
    def vcf_report(self):
        logging.info("Creating VCF")
        alt_seq = get_alt_seq(self.dir, self.out_nn, self.refpath)
        create_vcf(self.dir, self.thres, self.out_nn, self.refpath, self.rpath, self.rname, self.mapcmd, self.contig,
                   self.homo_t, self.het_t, self.minlen, self.depth, self.index2te, self.nv_cmd, alt_seq)
        logging.info("Creating HTML report")
        create_report(self.dir, self.contig, self.thres, self.rpath, self.refpath, self.rlendict, self.rname,
                      self.num_limit, self.ratio_limit)

    # Write intermediate data to file if debug mode
    def write2file(self, add_out):
        writer(os.path.join(self.dir, 'parse2.tsv'), add_out, self.debug)
        writer(os.path.join(self.dir, 'cluster.tsv'), self.out_nn, self.debug)

    # FASTA extractor
    def fasta_extract(self, qnames):
        fasta = os.path.join(self.dir, 'temp1.fa')
        out = open(os.path.join(self.dir, 'temp2.fa'), 'a')
        with open(fasta) as f:
            for line in f:
                try:
                    if qnames[line[1:].strip('\n')]:
                        line1 = line
                        line2 = next(f)
                        out.write(line1 + line2)
                except KeyError:
                    _ = next(f)
            out.close()
        # os.remove(fasta)


# Write to files
def writer(file_path, data_list, switch):
    if switch is True:
        w = open(file_path, 'w')
        w.write('\n'.join(data_list) + '\n')
        w.close()


# Create genome gap dictionary
def makegapdict(gap_path, contig_omit):
    gapdict = {}
    if gap_path is not None:
        rgapdata = open(gap_path, 'r').read().splitlines()
        for line in rgapdata:
            try:
                gapdict[line.split('\t')[0]].append([int(line.split('\t')[1]), int(line.split('\t')[2])])
            except KeyError:
                gapdict[line.split('\t')[0]] = [[int(line.split('\t')[1]), int(line.split('\t')[2])]]
    for contig in contig_omit:
        try:
            gapdict[contig].append(contig_omit[contig])
        except KeyError:
            gapdict[contig] = [contig_omit[contig]]
    if gapdict:  # Gap file is used or reference genome contigs contain invalid symbols
        logging.info('Gap dictionary successfully loaded')
        return gapdict
    else:
        logging.info('Gap dictionary not loaded.')
        return None

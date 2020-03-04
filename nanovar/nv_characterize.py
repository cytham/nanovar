"""
Method for SV characterization

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


import os
import logging
import random
from nanovar.nv_bam_parser import bam_parse
from nanovar.nv_detect_algo import sv_detect
from nanovar.nv_parser import entry_parser, align_info, breakpoint_parser
from nanovar.nv_cluster import sv_cluster
from nanovar.nv_nn import inference, svread_ovl
from nanovar.nv_cov_upper import ovl_upper
from nanovar.nv_vcf import create_vcf
from nanovar.nv_report import create_report


class VariantDetect:

    def __init__(self, wk_dir, bam, splitpct, minalign, filter_path, minlen, buff, model_path, total_gsize,
                 contig_len_dict, thres, read_path, read_name, ref_path, ref_name, map_cmd, mincov, homo_t, het_t, debug):
        self.dir = wk_dir
        self.bam = bam
        self.splitpct = splitpct
        self.minalign = minalign
        self.filter = filter_path
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
        self.basecov, self.maxovl, self.depth = 0, 0, 0
        self.total_out, self.total_subdata, self.out_nn, self.ins_out, self.out_rest = [], [], [], [], []
        self.rlendict, self.parse_dict = {}, {}
        # HTML SV table entry limit
        self.num_limit = 1000
        # HTML SV table SV ratio limit
        self.ratio_limit = 1
        self.maps = 0

    def bam_parse_detect(self):
        random.seed(1)
        self.total_subdata, self.total_out, self.basecov, self.parse_dict, self.rlendict, self.maps = bam_parse(self.bam,
                                                                                                                self.minlen,
                                                                                                                self.splitpct,
                                                                                                                self.minalign,
                                                                                                                self.dir,
                                                                                                                self.filter)

    def coverage_stats(self):
        # Obtaining upper cov limit and depth of coverage
        self.maxovl, self.depth = ovl_upper(self.gsize, self.contig, self.basecov, self.total_subdata, self.dir)
        # Report statistics of coverage
        logging.info("Genome size: %s bases" % str(self.gsize))
        logging.info("Mapped bases: %s bases" % str(self.basecov))
        logging.info("Depth of coverage: %sx" % str(self.depth))
        if self.depth < 4:
            logging.warning("Sequencing depth is less than 4x, output may not be comprehensive")
        logging.info("Read overlap upper limit: %s\n" % str(self.maxovl))

    def cluster_extract(self):
        logging.info("Clustering SV breakends")
        cluster_out = sv_cluster(self.total_subdata, self.total_out, self.buff, self.maxovl, self.mincov, self.contig, False)
        logging.info("Filtering INS and INV SVs")
        total_qnames = []
        for line in cluster_out:
            svtype = line.split('\t')[3].split(' ')[0]
            # cov = int(line.split('\t')[10])
            if svtype in ['S-Nov_Ins_bp', 'E-Nov_Ins_bp', 'Nov_Ins']:  # Some INS are actually DUP
                qnames = [line.split('\t')[0]]
                self.ins_out.append(self.parse_dict[line.split('\t')[8]])
                if line.split('\t')[11] != '.':
                    for mate in line.split('\t')[11].split(','):
                        qnames.append(mate[:-6])
                        self.ins_out.append(self.parse_dict[mate])
                total_qnames.extend(qnames)
            elif svtype in ['Inv', 'Inv(1)', 'Inv(2)']:  # For more precise INV bps
                qnames = [line.split('\t')[0]]
                if line.split('\t')[11] != '.':
                    for mate in line.split('\t')[11].split(','):
                        qnames.append(mate[:-6])
                total_qnames.extend(qnames)
            elif svtype not in ['S-Nov_Ins_bp', 'E-Nov_Ins_bp', 'Nov_Ins', 'Inv', 'Inv(1)', 'Inv(2)']:
                self.out_rest.append(line)
        total_qnames = set(total_qnames)
        self.fasta_extract(total_qnames)

    def cluster_nn(self, add_out):
        logging.info("Re-clustering INS/INV SVs and merging")
        # Merge old ins with new from hsblastn
        sub_out = self.ins_out + add_out
        cluster_out_ins = sv_cluster(self.total_subdata, sub_out, self.buff, self.maxovl, self.mincov, self.contig, True)
        new_cluster_out = self.out_rest + cluster_out_ins
        logging.info("Neural network inference")
        new_total_out = self.total_out + add_out
        if new_cluster_out:
            self.out_nn = inference(new_cluster_out, new_total_out, self.model)
        else:
            self.out_nn = []
        # Generate sv_overlap file
        svread_ovl(self.dir, self.out_nn)

    def parse_detect_hsb(self):
        random.seed(1)
        # Make gap dictionary
        gapdict = makegapdict(self.filter)
        temp1 = []
        chromocollect = []
        ovlt = 0.9
        sig_index = [0, 2, 4]
        data = open(self.bam, 'r').read().splitlines()
        data.append('null\tnull\tnull\tnull\tnull\tnull')
        nlines = len(data) - 1
        co = 0
        for i in range(nlines):
            if data[i].split('\t')[0] == data[i + 1].split('\t')[0]:  # Grouping alignments by read name
                temp1.append(align_info(data[i], self.rlendict))
                if data[i].split('\t')[1].strip() not in chromocollect:
                    chromocollect.append(data[i].split('\t')[1].strip())
            else:
                temp1.append(align_info(data[i], self.rlendict))
                if data[i].split('\t')[1].strip() not in chromocollect:
                    chromocollect.append(data[i].split('\t')[1].strip())
                co += 1
                # Parse entries and correct overlap alignments
                subdata = entry_parser(temp1, chromocollect, ovlt)
                for entry in subdata:
                    self.total_subdata.append('\t'.join(entry.split('\t')[0:5]))
                    # Add to base coverage
                    self.basecov += int(entry.split('\t')[2])
                # SV detection
                out1, out2 = sv_detect(subdata, self.splitpct, self.minalign, gapdict)
                if out2 == '':
                    pass
                else:
                    # Parse breakpoints
                    final = breakpoint_parser(out2, self.minlen, sig_index, co)
                    self.total_out.extend(final)
                temp1 = []
                chromocollect = []

    def vcf_report(self):
        logging.info("Creating VCF")
        create_vcf(self.dir, self.thres, self.out_nn, self.refpath, self.rpath, self.rname, self.mapcmd, self.contig,
                   self.homo_t, self.het_t)
        logging.info("Creating HTML report")
        create_report(self.dir, self.contig, self.thres, self.rpath, self.refpath, self.rlendict, self.rname,
                      self.num_limit, self.ratio_limit)

    # Write intermediate data to file if debug mode
    def write2file(self, add_out):
        writer(os.path.join(self.dir, 'subdata.tsv'), self.total_subdata, self.debug)
        writer(os.path.join(self.dir, 'parse1.tsv'), self.total_out, self.debug)
        writer(os.path.join(self.dir, 'parse2.tsv'), add_out, self.debug)
        writer(os.path.join(self.dir, 'cluster.tsv'), self.out_nn, self.debug)

    # Parse data and detect SVs
    def parse_detect(self, total_lines, contig_collect, seed, gapdict, ovlt, sig_index):
        lines_sort = sorted(sorted(total_lines, key=lambda x: x[1], reverse=True), key=lambda y: y[0])
        temp1 = [tup[2] for tup in lines_sort]
        # Parse entries and correct overlap alignments
        subdata = entry_parser(temp1, contig_collect, ovlt)
        for entry in subdata:
            self.total_subdata.append('\t'.join(entry.split('\t')[0:5]))
            # Add to base coverage
            self.basecov += int(entry.split('\t')[2])
        # SV detection
        out1, out2 = sv_detect(subdata, self.splitpct, self.minalign, gapdict)
        if out1 == '' and out2 == '':
            pass
        else:
            # Parse breakpoints
            final = breakpoint_parser(out2, self.minlen, sig_index, seed)
            self.total_out.extend(final)
            for i in final:
                self.parse_dict[i.split('\t')[8]] = i

    # FASTA extractor
    def fasta_extract(self, qnames):
        fasta = os.path.join(self.dir, 'temp1.fa')
        out = open(os.path.join(self.dir, 'temp2.fa'), 'a')
        with open(fasta) as f:
            for line in f:
                if line.strip('\n') in qnames:
                    line1 = line
                    line2 = next(f)
                    out.write('>' + line1 + line2)
                else:
                    line2 = next(f)
            out.close()
        os.remove(fasta)


# Write to files
def writer(file_path, data_list, switch):
    if switch is True:
        w = open(file_path, 'w')
        w.write('\n'.join(data_list) + '\n')
        w.close()


# Create genome gap dictionary
def makegapdict(gap_path):
    gapdict = {}
    if gap_path is not None:
        rgapdata = open(gap_path, 'r').read().splitlines()
        n = len(rgapdata)
        for i in range(n):  # Define gap dict list
            gapdict[rgapdata[i].split('\t')[0]] = []
        for i in range(n):  # Making gap dictionary
            gapdict[rgapdata[i].split('\t')[0]].append([int(rgapdata[i].split('\t')[1]), int(rgapdata[i].split('\t')[2])])
        logging.info('Gap dictionary successfully loaded')
        return gapdict
    else:
        logging.info('Gap dictionary not loaded.')
        return None

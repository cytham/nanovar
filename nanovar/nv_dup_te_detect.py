"""
Duplication and transposable element analyzer. 

Copyright (C) 2023 Tham Cheng Yong

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
import pysam
import logging
from pybedtools import BedTool
# from subprocess import Popen, PIPE, STDOUT
from nanovar.nv_align import align_mm

def dup_te_analyzer(wk_dir, out_nn, total_out1, score_threshold, ref_dir, ref_path, mm, threads, data_type, st, debug):
    index2reads = get_index2reads(out_nn, score_threshold)
    region_dict, ins_coord_dict, dups = get_region_dict(total_out1)
    read2index = get_ins_sequence(index2reads, region_dict, wk_dir)
    te_ref_file = os.path.join(ref_dir, 'hg38_L1_Alu.fa')
    reads_file = os.path.join(wk_dir, 'ins_seq.fa')
    _, bam_te = align_mm(te_ref_file, reads_file, wk_dir, 'ins_seq', 'te', threads, mm, data_type, st)
    index2te = parse_bam_te(bam_te, wk_dir, read2index)
    unmap_reads_file = os.path.join(wk_dir, 'ins_seq_te-unmap.fa')
    _, bam_ref = align_mm(ref_path, unmap_reads_file, wk_dir, 'unmap', 'ref', threads, mm, data_type, st)
    index2dup, dup_coord = parse_bam_ref(bam_ref, ins_coord_dict, read2index)
    index2dup = update_dups(dups, index2dup, index2reads)
    dup_indexes, _ = filter_dup(index2dup, index2reads, threshold=0.75)
    new_out_nn = update_out_nn(dup_indexes, out_nn)
    if not debug:
        os.remove(reads_file)
        os.remove(unmap_reads_file)
        os.remove(bam_te)
        os.remove(bam_ref)
        
    return index2te, new_out_nn

def get_index2reads(out_nn, score_threshold):
    # Convert score threshold to nn probability
    nn_score = 1-10**(score_threshold/-10)
    index2reads = {}
    for i in out_nn:
        if i.split('\t')[3].split(' ')[0] in ['Nov_Ins', 'S-Nov_Ins_bp', 'E-Nov_Ins_bp']:
            if float(i.split('\t')[13]) >= nn_score:
                ind = i.split('\t')[6].split('~')[0]
                reads = i.split('\t')[8] + ',' + i.split('\t')[11]
                index2reads[ind] = reads
    return index2reads

def get_region_dict(total_out):
    region_dict = {}
    ins_coord_dict = {}
    dups = {}
    for i in total_out:
        sv_type = i.split('\t')[3]
        if sv_type[:3] in ('S-N', 'E-N', 'Nov'):
            readindex = i.split('\t')[8]
            query_ins_start = int(i.split('\t')[1])
            ins_size = float(i.split('\t')[3].split(' ')[1].split('~')[0])
            chrm = i.split('\t')[6].split('~')[1].split(':')[0]
            coord = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            ins_coord_dict[readindex] = [chrm, coord-ins_size, coord+ins_size]
            if sv_type[:3] == 'S-N':
                region_dict[readindex] = (1, ins_size)
            elif sv_type[:3] == 'E-N':
                region_dict[readindex] = (query_ins_start, query_ins_start + ins_size)
            elif sv_type[:3] == 'Nov':
                region_dict[readindex] = (query_ins_start, query_ins_start + ins_size)
        elif sv_type[:5] == 'TDupl':
            readindex = i.split('\t')[8]
            dups[readindex] = ''
    return region_dict, ins_coord_dict, dups

def get_ins_sequence(index2reads, region_dict, wk_dir):
    bed_str = ''
    read2index = {}
    fai_path = os.path.join(wk_dir, 'temp1.fa.fai')
    fai_names = {i.split('\t')[0]: '' for i in open(fai_path, 'r').read().splitlines()}
    missed_reads = 0
    for i in index2reads:
        for r in index2reads[i].split(','):
            if r.split('~')[0] in fai_names:
                try:
                    left = str(region_dict[r.strip('\n')][0])
                    right = str(int(region_dict[r.strip('\n')][1]))
                    if r.split('~')[0] + ':' + left + '-' + right not in read2index:
                        bed_str += r.split('~')[0] + '\t' + left + '\t' + right + '\t' + r.strip('\n') + '\n'
                        read2index[r.strip('\n')] = i
                except KeyError:  ## Some reads might not be previously assigned as INS (e.g. TDupl), hence not picked up by region_dict
                    pass
            else:
                missed_reads += 1
    log_missing_reads(missed_reads)
    bed = BedTool(bed_str, from_string=True)
    # Extract fasta for each insertion region
    fasta_path = os.path.join(wk_dir, 'temp1.fa')
    fasta = bed.sequence(fi=fasta_path, name=True)
    _ = fasta.save_seqs(os.path.join(wk_dir, 'ins_seq.fa'))
    return read2index

def log_missing_reads(missed_reads):
    if missed_reads > 0:
        logging.warning("Could not analyze INS sequence of %i reads due to absence of FASTA. Please check if all reads in BAM have a primary alignment." % missed_reads)

def parse_bam_te(bam, wk_dir, read2index, L1_90_len=5500, Alu_90_len=270):
    save = pysam.set_verbosity(0)
    sam = pysam.AlignmentFile(bam, "rb")
    pysam.set_verbosity(save)
    index2te = {}
    unmap_fasta = open(os.path.join(wk_dir, 'ins_seq_te-unmap.fa'), 'w')
    for seg in sam:
        flag = seg.flag
        qname = seg.query_name.split('::')[0]
        if flag == 4:
            unmap_fasta.write('>' + seg.query_name + '\n' + seg.query_sequence + '\n')
        elif flag in (0, 16):
            ref_aligned_len = seg.reference_end - seg.reference_start
            if seg.reference_name.rsplit('_', 1)[1].startswith('L1'):
                if int(ref_aligned_len) >= L1_90_len:  # 90% of L1
                    pass
            elif seg.reference_name.rsplit('_', 1)[1].startswith('Alu'):
                if int(ref_aligned_len) >= Alu_90_len:  # 90% of Alu
                    pass
            else:
                continue
            if read2index[qname] in index2te:
                index2te[read2index[qname]].add(seg.reference_name.rsplit('_', 1)[1])
            else:
                index2te[read2index[qname]] = set()
                index2te[read2index[qname]].add(seg.reference_name.rsplit('_', 1)[1])
    unmap_fasta.close()
    return index2te

def parse_bam_ref(bam, ins_coord_dict, read2index):
    save = pysam.set_verbosity(0)
    sam = pysam.AlignmentFile(bam, "rb")
    pysam.set_verbosity(save)
    index2dup = {}
    dup_coord = {}
    for seg in sam:
        flag = seg.flag
        if flag in (0, 16):
            qname = seg.query_name.split('::')[0]
            dup, coord = overlap_dup(ins_coord_dict[qname], seg)
            if dup:
                dup_coord[qname] = coord
                if read2index[qname] in index2dup:
                    index2dup[read2index[qname]] += 1
                else:
                    index2dup[read2index[qname]] = 1
    return index2dup, dup_coord

def overlap_dup(supposed_coord, seg, buffer=20):
    if seg.reference_name == supposed_coord[0]:
        if supposed_coord[1] - buffer <= seg.reference_start and seg.reference_end <= supposed_coord[2] + buffer:
            return True, seg.reference_name + ':' + str(seg.reference_start) + '-' +  str(seg.reference_end)
        else:
            return False, ''
    else:
        return False, ''

# Add originally TDupl reads
def update_dups(dups, index2dup, index2reads):
    for i in index2reads:
        for r in index2reads[i].split(','):
            if r in dups:
                try:
                    index2dup[i] += 1
                except KeyError:
                    pass
    return index2dup

# Find TDupls based on proportions of dup reads
def filter_dup(index2dup, index2reads, threshold=0.75):
    dup_indexes = {}
    fail_dups = {}
    for i in index2dup:
        if index2dup[i]/len(index2reads[i].split(',')) >= threshold:
            dup_indexes[i] = round(index2dup[i]/len(index2reads[i].split(',')), 3)
        else:
            # fail_dups[i] = round(index2dup[i]/len(index2reads[i].split(',')), 3)
            fail_dups[i] = ''
    return dup_indexes, fail_dups

def update_out_nn(dup_indexes, out_nn, dup_coord):
    new_out_nn = []
    for i in out_nn:
        if i.split('\t')[6].split('~')[0] in dup_indexes:
            index_coord = i.split('\t')[6].split('~')[0] + '~' + dup_coord[i.split('\t')[8]]
            new_out_nn.append('\t'.join(i.split('\t')[0:3]) + '\tTDupl 99.99~' + i.split('\t')[3].split('~')[1] + '\t' + '\t'.join(i.split('\t')[4:6]) + '\t' + index_coord + '\t' + '\t'.join(i.split('\t')[7:]))
        else:
            new_out_nn.append(i)
    return new_out_nn

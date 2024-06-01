"""
Extract alternative sequences

Copyright (C) 2024 Tham Cheng Yong

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
from pybedtools import BedTool
from Bio.Seq import Seq

def get_alt_seq(wk_dir, out_nn, ref_path):
    bed_str = ''
    ins_id_seq = {}
    ins_read_seq = get_ins_seq(wk_dir)
    for line in out_nn:
       sv_type = line.split('\t')[3].split(' ')[0]
       sv_id = line.split('\t')[6].split('~')[0]
       strand = line.split('\t')[3].split('~')[1].split(',')[0]
       bed_line = make_bed(line, sv_type, sv_id)
       bed_str += bed_line
       if sv_type in ['Nov_Ins', 'E-Nov_Ins_bp', 'S-Nov_Ins_bp']:
           ins_seq = ins_seq_except(ins_read_seq, line.split('\t')[8])
           if strand == '+':
               ins_id_seq[sv_id] = ins_seq
           elif strand == '-':
               ins_id_seq[sv_id] = str(Seq(ins_seq).reverse_complement())
           else:
               raise Exception('{} strand symbol invalid'.format(strand))
    bed = BedTool(bed_str, from_string=True)
    fasta = bed.sequence(fi=ref_path, nameOnly=True)
    alt_seq = {}
    with open(fasta.seqfn) as f:
        for r in f:
            alt_seq[r.strip('>\n')] = next(f).upper().strip()
    for i in ins_id_seq:
        alt_seq[i] = alt_seq[i] + ins_id_seq[i]
    return alt_seq

def make_bed(line, sv_type, sv_id):
    data = line.split('\t')
    chrm = data[6].split('~')[1].split(':')[0]
    if sv_type in ['Inter-Ins(1)', 'Inter-Ins(2)', 'InterTx']:
        end = int(data[6].split('~')[1].split(':')[1])
        start = end - 1
    elif sv_type == 'Del':
        start = int(data[6].split('~')[1].split(':')[1].split('-')[0]) - 1
        end = data[6].split('~')[1].split(':')[1].split('-')[1]
    else:  # ['Nov_Ins', 'E-Nov_Ins_bp', 'S-Nov_Ins_bp', 'Inv', 'Inv(1)', 'Inv(2)', 'TDupl', 'Intra-Ins', 'Intra-Ins(1)', 'Intra-Ins(2)']
        end = int(data[6].split('~')[1].split(':')[1].split('-')[0])
        start = end - 1
    start = min(start, 0)
    bed_line = '\t'.join([chrm, str(start), str(end), sv_id]) + '\n'
    return bed_line

def get_ins_seq(wk_dir):
    with open(os.path.join(wk_dir, 'ins_seq.fa')) as f:
        ins_read_seq = {}
        for r in f:
            id = r.split('::')[0].strip('>')
            ins_read_seq[id] = next(f).strip()
        return ins_read_seq

def ins_seq_except(ins_read_seq, read_id):
    try:
        return ins_read_seq[read_id]
    except KeyError:  # Not found in ins_seq.fa, below score threshold
        return 'N'


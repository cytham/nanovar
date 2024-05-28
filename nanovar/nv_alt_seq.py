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

def get_alt_seq(out_nn, ref_path):
    bed_str = ''
    for line in out_nn:
       bed_line = make_bed(line)
       bed_str += bed_line
    bed = BedTool(bed_str, from_string=True)
    fasta = bed.sequence(fi=ref_path, nameOnly=True)
    alt_seq = {}
    with open(fasta.seqfn) as f:
        for r in f:
            alt_seq[r.strip('>\n')] = next(f).upper().strip()
    return alt_seq

def make_bed(line):
    data = line.split('\t')
    sv_type = data[3].split(' ')[0]
    sv_id = data[6].split('~')[0]
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
    bed_line = '\t'.join([chrm, str(start), str(end), sv_id]) + '\n'
    return bed_line

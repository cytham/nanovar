"""
Transposable element analyzer.

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
from pybedtools import BedTool
from subprocess import Popen, PIPE, STDOUT


def te_analyzer(wk_dir, out_nn, total_out1, total_out2, score_threshold, ref_dir, hsb, mdb, threads, debug):
    # Convert score threshold to nn probability
    nn_score = 1-10**(score_threshold/-10)
    index2reads = {}
    for i in out_nn:
        if i.split('\t')[3].split(' ')[0] in ['Nov_Ins', 'S-Nov_Ins_bp', 'E-Nov_Ins_bp']:
            if float(i.split('\t')[13]) >= nn_score:
                ind = i.split('\t')[6].split('~')[0]
                reads = i.split('\t')[8] + ',' + i.split('\t')[11]
                index2reads[ind] = reads
    region_dict = {}
    for i in total_out1:
        if i.split('\t')[3][:3] == 'S-N':
            region_dict[i.split('\t')[8]] = (1, float(i.split('\t')[3].split(' ')[1].split('~')[0]))
        elif i.split('\t')[3][:3] == 'E-N':
            region_dict[i.split('\t')[8]] = (int(i.split('\t')[1]), int(i.split('\t')[1]) +
                                             float(i.split('\t')[3].split(' ')[1].split('~')[0]))
        elif i.split('\t')[3][:7] == 'Nov_Ins':
            region_dict[i.split('\t')[8]] = (int(i.split('\t')[1]), int(i.split('\t')[1]) +
                                             float(i.split('\t')[3].split(' ')[1].split('~')[0]))
        else:
            continue
    for i in total_out2:
        if i.split('\t')[3][:3] == 'S-N':
            region_dict[i.split('\t')[8]] = (1, float(i.split('\t')[3].split(' ')[1].split('~')[0]))
        elif i.split('\t')[3][:3] == 'E-N':
            region_dict[i.split('\t')[8]] = (int(i.split('\t')[1]), int(i.split('\t')[1]) +
                                             float(i.split('\t')[3].split(' ')[1].split('~')[0]))
        elif i.split('\t')[3][:7] == 'Nov_Ins':
            region_dict[i.split('\t')[8]] = (int(i.split('\t')[1]), int(i.split('\t')[1]) +
                                             float(i.split('\t')[3].split(' ')[1].split('~')[0]))
        else:
            continue
    bed_str = ''
    read2index = {}
    for i in index2reads:
        for r in index2reads[i].split(','):
            try:
                left = str(region_dict[r.strip('\n')][0])
                right = str(int(region_dict[r.strip('\n')][1]))
                if r.split('~')[0] + ':' + left + '-' + right not in read2index:
                    bed_str += r.split('~')[0] + '\t' + left + '\t' + right + '\n'
                    read2index[r.split('~')[0] + ':' + left + '-' + right] = i
            except KeyError:
                pass
    bed = BedTool(bed_str, from_string=True)
    # Extract fasta for each insertion region
    fasta_path = os.path.join(wk_dir, 'temp1.fa')
    fasta = bed.sequence(fi=fasta_path)
    _ = fasta.save_seqs(os.path.join(wk_dir, 'ins_seq.fa'))
    # HS-BLASTN alignment
    ref = os.path.join(ref_dir, 'hg38_L1_Alu.fa')
    # Index TE reference (makeblastdb)
    process = Popen([mdb, '-in', ref, '-input_type', 'fasta', '-dbtype', 'nucl'], universal_newlines=True,
                    stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: makeblastdb failed (2)")
        raise Exception("Error: makeblastdb failed (2), see log")
    # Index TE reference (hs-blastn)
    process = Popen([hsb, 'index', ref], universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: hs-blastn index failed (2)")
        raise Exception("Error: hs-blastn index failed (2), see log")
    obinary_path = os.path.join(ref_dir, 'hg38_L1_Alu.counts.obinary')
    read_path = os.path.join(wk_dir, 'ins_seq.fa')
    out_path = os.path.join(wk_dir, 'temp-te-blast.tab')
    process = Popen([hsb + ' align -db ' + ref + ' -window_masker_db ' + obinary_path + ' -query ' + read_path + ' -out '
                     + out_path + ' -outfmt 6 -num_threads ' + str(threads) + ' -max_target_seqs 3 -gapopen 0 '
                    '-gapextend 4 -penalty -3 -reward 2'],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: hs-blastn alignment failed (2)")
        raise Exception("Error: hs-blastn alignment failed (2), see log")
    index2te = {}
    with open(out_path) as f:
        for line in f:
            line = line.split('\t')
            if line[1].rsplit('_', 1)[1].startswith('L1'):
                if int(line[3]) >= 5500:  # 90% of L1
                    if read2index[line[0]] in index2te:
                        index2te[read2index[line[0]]].add(line[1].rsplit('_', 1)[1])
                    else:
                        index2te[read2index[line[0]]] = set()
                        index2te[read2index[line[0]]].add(line[1].rsplit('_', 1)[1])
            elif line[1].rsplit('_', 1)[1].startswith('Alu'):
                if int(line[3]) >= 270:  # 90% of Alu
                    if read2index[line[0]] in index2te:
                        index2te[read2index[line[0]]].add(line[1].rsplit('_', 1)[1])
                    else:
                        index2te[read2index[line[0]]] = set()
                        index2te[read2index[line[0]]].add(line[1].rsplit('_', 1)[1])
    if not debug:
        os.remove(read_path)
        os.remove(out_path)
        
    return index2te


def log_subprocess(out):
    for line in iter(out.readline, ''):
        if line != '\n':
            logging.debug(line.strip())

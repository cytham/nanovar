"""
Parsing algorithm.

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

import random
import re
import ast


# Parse alignment entries
def entry_parser(temp1, chromocollect, overlap_tol):
    temp2 = []
    num_chr = len(chromocollect)  # Total number of chromosomes aligned
    num_align = len(temp1)  # Total number of alignments
    while len(temp1) != 0:
        # Lock-in lead alignment with the highest bitscore
        temp2.append(temp1[0] + '\tn=' + str(num_align) + '\t' + str(num_chr) + 'chr')
        # leadrange = [int(temp1[0].split('\t')[5]), int(temp1[0].split('\t')[5]) + int(temp1[0].split('\t')[6])]
        leadstart = int(temp1[0].split('\t')[5])
        leadend = int(temp1[0].split('\t')[5]) + int(temp1[0].split('\t')[6])
        del temp1[0]  # Remove lead alignment entry
        j = len(temp1)  # Recalculate number of remaining entries
        temp3 = []
        for p in range(j):
            qstart = int(temp1[p].split('\t')[5])
            qend = int(temp1[p].split('\t')[5]) + int(temp1[p].split('\t')[6])
            queryrange = range(qstart, qend + 1)
            leadintersect = range(max(leadstart, qstart), min(leadend, qend) + 1)
            if len(leadintersect) == 0:  # No intersection with lead
                temp3.append(temp1[p])
            else:
                qulen = len(queryrange)
                leadintlen = len(leadintersect)
                if float(leadintlen) / qulen > overlap_tol or (qulen - leadintlen) <= 50:
                    # If overlap len is more than 90% of query alignment length or
                    # after-trimmed length <=50 bp, lead wins and query is omitted. Otherwise, do trimming.
                    continue  # Entry omitted
                else:  # If overlap is tolerated, then do trimming
                    new = ''
                    strand = str(temp1[p].split('\t')[7])
                    # Alter query and subject ranges and adjust bitscore as a proportion of alignment length
                    if qstart == min(leadintersect):  # Left overlap
                        if strand == '+':
                            new = temp1[p].split('\t')[0] + '\t' + str(
                                int(temp1[p].split('\t')[1]) + leadintlen) + '\t' + str(
                                int(temp1[p].split('\t')[2]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[3:5]) + '\t' + str(
                                int(temp1[p].split('\t')[5]) + leadintlen) + '\t' + str(
                                int(temp1[p].split('\t')[6]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[7:10]) + '\t' + str(
                                round(float(temp1[p].split('\t')[10]) * (float(qulen - leadintlen) / qulen),
                                      2)) + '\t' + '\t'.join(temp1[p].split('\t')[11:])
                        elif strand == '-':
                            new = temp1[p].split('\t')[0] + '\t' + temp1[p].split('\t')[1] + '\t' + str(
                                int(temp1[p].split('\t')[2]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[3:5]) + '\t' + str(
                                int(temp1[p].split('\t')[5]) + leadintlen) + '\t' + str(
                                int(temp1[p].split('\t')[6]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[7:10]) + '\t' + str(
                                round(float(temp1[p].split('\t')[10]) * (float(qulen - leadintlen) / qulen),
                                      2)) + '\t' + '\t'.join(temp1[p].split('\t')[11:])
                    elif qend == max(leadintersect):  # Right overlap
                        if strand == '+':
                            new = temp1[p].split('\t')[0] + '\t' + temp1[p].split('\t')[1] + '\t' + str(
                                int(temp1[p].split('\t')[2]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[3:5]) + '\t' + temp1[p].split('\t')[5] + '\t' + str(
                                int(temp1[p].split('\t')[6]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[7:10]) + '\t' + str(
                                round(float(temp1[p].split('\t')[10]) * (float(qulen - leadintlen) / qulen),
                                      2)) + '\t' + '\t'.join(temp1[p].split('\t')[11:])
                        elif strand == '-':
                            new = temp1[p].split('\t')[0] + '\t' + str(
                                int(temp1[p].split('\t')[1]) + leadintlen) + '\t' + str(
                                int(temp1[p].split('\t')[2]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[3:5]) + '\t' + temp1[p].split('\t')[5] + '\t' + str(
                                int(temp1[p].split('\t')[6]) - leadintlen) + '\t' + '\t'.join(
                                temp1[p].split('\t')[7:10]) + '\t' + str(
                                round(float(temp1[p].split('\t')[10]) * (float(qulen - leadintlen) / qulen),
                                      2)) + '\t' + '\t'.join(temp1[p].split('\t')[11:])
                    else:  # If lead lies entirely within intersection, omit entry
                        continue
                    if new != '' and int(new.split('\t')[2]) > 0:
                        temp3.append(new)
        temp1 = temp3
    sortdict = {}
    for k in temp2:
        sortdict[int(k.split('\t')[5])] = k
    return [value for (key, value) in sorted(sortdict.items())]


# Extract alignment info
def align_info(line, rlendict):
    qid = line.split('\t')[0]
    sid = line.split('\t')[1]
    piden = line.split('\t')[2]
    nmismatch = line.split('\t')[4]
    ngapopen = line.split('\t')[5]
    qstart = line.split('\t')[6]
    qstretch = int(line.split('\t')[7]) - int(line.split('\t')[6])
    sstart = min(int(line.split('\t')[8]), int(line.split('\t')[9]))
    sstretch = abs(int(line.split('\t')[9]) - int(line.split('\t')[8]))
    evalue = line.split('\t')[10]
    bitscore = line.split('\t')[11].strip('\n')
    if line.split('\t')[8] <= line.split('\t')[9]:
        strand = "+"
    else:
        strand = "-"
    rlen = rlendict[qid]
    entry = sid + '\t' + str(sstart) + '\t' + str(sstretch) + '\t+\t' + qid + '\t' + qstart + '\t' + str(qstretch) + '\t' + \
        strand + '\t' + str(rlen) + '\t' + evalue + '\t' + bitscore + '\t' + piden + '\t' + nmismatch + '\t' + ngapopen
    return entry


# Parse SV breakpoints
def breakpoint_parser(out, minlen, sig_index, seed):
    final = []
    ran = "01234567890ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    ran_len = 5
    read_name = out.split('\t')[0]
    random.seed(seed)
    if out.split('\t')[17] == '':
        complex_sv = 0
    else:
        complex_sv = len(out.split('\t')[17].split(','))
    # if int(out.split('\t')[17].count(',')) >= 1:
    #     multi_sv = len(out.split('\t')[18].split(','))
    # else:
    #     multi_sv = 0
    del_count = del_counter = int(out.split('\t')[19].count('Del'))
    iso_ins_count = iso_ins_counter = int(out.split('\t')[19].count('S-Nov_Ins')) + int(out.split('\t')[19].count('E-Nov_Ins'))
    ins_count = ins_counter = int(out.split('\t')[19].count('Nov_Ins'))
    tdup_count = tdup_counter = int(out.split('\t')[19].count('TDupl'))
    inv_count = inv_counter = int(out.split('\t')[19].count('Inv'))
    intra_ins_count = intra_ins_counter = int(out.split('\t')[19].count('Intra-Ins'))
    inter_ins_count = inter_ins_counter = int(out.split('\t')[19].count('Inter-Ins'))
    inter_tx_count = inter_tx_counter = int(out.split('\t')[19].count('InterTx'))
    strandness = out.split('\t')[15].split(',')
    bps = out.split('\t')[19].split(',')
    bps_no = len(bps)
    querymap = out.split('\t')[20]
    nchrom = out.split('\t')[1]
    evaluelist = out.split('\t')[4].split(',')
    bitscorelist = out.split('\t')[5].split(',')
    short = out.split('\t')[16]
    piden = out.split('\t')[21].split(',')
    mismatch = out.split('\t')[22].split(',')
    gap_ratio = out.split('\t')[23].split(',')
    sign = out.split('\t')[3]
    realnmaps = int(out.split('\t')[2].split(' ')[0])
    # if realnmaps >= 5:
    #     final = []
    #     return final
    for j in range(bps_no):
        bp_name = bps[j].split(':')[0]
        bp_range = re.sub('^.*:', '', str(bps[j]))
        bp1 = int(bp_range.split('-')[0])
        bp2 = int(bp_range.split('-')[1])
        signature, strands = getsignature(bp1, bp_name, querymap, sign, realnmaps, evaluelist, bitscorelist, bps_no,
                                          complex_sv, nchrom, short, piden, mismatch, gap_ratio, strandness, sig_index)
        if bp_name == 'S-Nov_Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[0]
            coord = out.split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[1]
            s_nov_ins_size = out.split('\t')[7].split(' ')[0].split(',')[int(iso_ins_count - iso_ins_counter)].split('ns')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord)
            uniqname = read_name + '~' + bp_uname
            if int(s_nov_ins_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp2 - 10) + '\t' + str(bp2) + '\t' + str(bp_name + '_bp') + ' ' +
                    str(s_nov_ins_size) + '~' + str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                )
            iso_ins_counter = iso_ins_counter - 1
        elif bp_name == 'E-Nov_Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[0]
            coord = out.split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[1]
            e_nov_ins_size = out.split('\t')[7].split(' ')[0].split(',')[int(iso_ins_count - iso_ins_counter)].split('ns')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord)
            uniqname = read_name + '~' + bp_uname
            if int(e_nov_ins_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp1 + 10) + '\t' + str(bp_name + '_bp') + ' ' +
                    str(e_nov_ins_size) + '~' + str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                )
            iso_ins_counter = iso_ins_counter - 1
        elif bp_name == 'Del':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[6].split(' ')[1].split(',')[int(del_count - del_counter)].split(':')[0]
            coord1 = int(out.split('\t')[6].split(' ')[1].split(',')[int(del_count - del_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[6].split(' ')[1].split(',')[int(del_count - del_counter)].split(':')[1].split('-')[1])
            del_size = out.split('\t')[6].split(' ')[0].split(',')[int(del_count - del_counter)].split('l')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1, coord2)) + '-' + str(max(coord1, coord2))
            uniqname = read_name + '~' + bp_uname
            if int(del_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(del_size) + '~' +
                    str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                    str(uniqname) + '\t' + str(signature)
                )
            del_counter = del_counter - 1
        elif bp_name == 'Nov_Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[8].split(' ')[1].split(',')[int(ins_count - ins_counter)].split(':')[0]
            coord1 = int(out.split('\t')[8].split(' ')[1].split(',')[int(ins_count - ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[8].split(' ')[1].split(',')[int(ins_count - ins_counter)].split(':')[1].split('-')[1])
            nov_ins_size = out.split('\t')[8].split(' ')[0].split(',')[int(ins_count - ins_counter)].split('ns')[1]
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1, coord2)) + '-' + str(max(coord1, coord2))
            uniqname = read_name + '~' + bp_uname
            if int(nov_ins_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp1 + 10) + '\t' + str(bp_name) + ' ' +
                    str(nov_ins_size) + '~' + str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str('Nov1') + '\t' + str(uniqname) + '\t' + str(signature)
                )
            ins_counter = ins_counter - 1
        elif bp_name == 'TDupl':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[10].split(' ')[1].split(',')[int(tdup_count - tdup_counter)].split(':')[0]
            coord1 = int(out.split('\t')[10].split(' ')[1].split(',')[int(tdup_count - tdup_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[10].split(' ')[1].split(',')[int(tdup_count - tdup_counter)].split(':')[1].split('-')[1])
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1, coord2)) + '-' + str(max(coord1, coord2))
            uniqname = read_name + '~' + bp_uname
            final.append(
                str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '~' +
                str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                str(uniqname) + '\t' + str(signature)
            )
            tdup_counter = tdup_counter - 1
        elif bp_name == 'Inv':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[0]
            coord1 = int(out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[1])
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1, coord2)) + '-' + str(max(coord1, coord2))
            uniqname = read_name + '~' + bp_uname
            final.append(
                str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '~' +
                str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                str(uniqname) + '\t' + str(signature)
            )
            inv_counter = inv_counter - 1
        elif bp_name == 'Inv(1)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[0]
            coord1 = int(out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[1])
            inv_size = out.split('\t')[11].split(' ')[0].split(',')[int(inv_count - inv_counter)].split('nv')[1]
            pair = 'Inv1'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1, coord2)) + '-' + str(max(coord1, coord2))
            uniqname = read_name + '~' + bp_uname
            if int(inv_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(inv_size) + '~' +
                    str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                    str(uniqname) + '\t' + str(signature)
                )
            inv_counter = inv_counter - 1
        elif bp_name == 'Inv(2)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[0]
            coord1 = int(out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[1])
            inv_size = out.split('\t')[11].split(' ')[0].split(',')[int(inv_count - inv_counter)].split('nv')[1]
            pair = 'Inv2'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1, coord2)) + '-' + str(max(coord1, coord2))
            uniqname = read_name + '~' + bp_uname
            if int(inv_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(inv_size) + '~' +
                    str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                    str(uniqname) + '\t' + str(signature)
                )
            inv_counter = inv_counter - 1
        elif bp_name == 'Intra-Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[0]
            coord1 = int(out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count -
                                                                          intra_ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count -
                                                                          intra_ins_counter)].split(':')[1].split('-')[1])
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord1) + '-' + str(coord2)
            uniqname = read_name + '~' + bp_uname
            final.append(
                str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '~' +
                str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                str(uniqname) + '\t' + str(signature)
            )
            intra_ins_counter = intra_ins_counter - 1
        elif bp_name == 'Intra-Ins(1)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[0]
            coord1 = int(out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count -
                                                                          intra_ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count -
                                                                          intra_ins_counter)].split(':')[1].split('-')[1])
            intrains_size = out.split('\t')[12].split(' ')[0].split(',')[int(intra_ins_count - intra_ins_counter)].split(')')[1]
            pair = 'Intra1'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord1) + '-' + str(coord2)
            uniqname = read_name + '~' + bp_uname
            if int(intrains_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' +
                    str(intrains_size) + '~' + str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                )
            intra_ins_counter = intra_ins_counter - 1
        elif bp_name == 'Intra-Ins(2)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[0]
            coord1 = int(out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count -
                                                                          intra_ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(out.split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count -
                                                                          intra_ins_counter)].split(':')[1].split('-')[1])
            intrains_size = out.split('\t')[12].split(' ')[0].split(',')[int(intra_ins_count - intra_ins_counter)].split(')')[1]
            pair = 'Intra2'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord1) + '-' + str(coord2)
            uniqname = read_name + '~' + bp_uname
            if int(intrains_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' +
                    str(intrains_size) + '~' + str(strands) + '\t' + str(chrom) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                )
            intra_ins_counter = intra_ins_counter - 1
        elif bp_name == 'Inter-Ins(1)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom1 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                      inter_ins_counter)].split('~')[0].split(':')[0]
            coord1_1 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                        inter_ins_counter)].split('~')[0].split(':')[1]
            chrom2 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                      inter_ins_counter)].split('~')[1].split(':')[0]
            coord2_1 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                        inter_ins_counter)].split('~')[1].split(':')[1]
            interins_size = out.split('\t')[13].split(' ')[0].split(',')[int(inter_ins_count - inter_ins_counter)].split(')')[1]
            pair = 'Inter1'
            uniqidx = bp_uname + '~' + str(chrom1) + ':' + str(coord1_1) + '~' + str(chrom2) + ':' + str(coord2_1)
            uniqname = read_name + '~' + bp_uname
            if int(interins_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' +
                    str(interins_size) + '~' + str(strands) + '\t' + str(chrom1) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                )
            inter_ins_counter = inter_ins_counter - 1
        elif bp_name == 'Inter-Ins(2)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom1 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                      inter_ins_counter)].split('~')[0].split(':')[0]
            coord1_1 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                        inter_ins_counter)].split('~')[0].split(':')[1]
            chrom2 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                      inter_ins_counter)].split('~')[1].split(':')[0]
            coord2_1 = out.split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count -
                                                                        inter_ins_counter)].split('~')[1].split(':')[1]
            interins_size = out.split('\t')[13].split(' ')[0].split(',')[int(inter_ins_count - inter_ins_counter)].split(')')[1]
            pair = 'Inter2'
            uniqidx = bp_uname + '~' + str(chrom1) + ':' + str(coord1_1) + '~' + str(chrom2) + ':' + str(coord2_1)
            uniqname = read_name + '~' + bp_uname
            if int(interins_size) >= minlen:
                final.append(
                    str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' +
                    str(interins_size) + '~' + str(strands) + '\t' + str(chrom1) + '\t' + str('-') + '\t' +
                    str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                )
            inter_ins_counter = inter_ins_counter - 1
        elif bp_name == 'InterTx':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom1 = out.split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count -
                                                                      inter_tx_counter)].split('~')[0].split(':')[0]
            coord1_1 = out.split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count -
                                                                        inter_tx_counter)].split('~')[0].split(':')[1]
            chrom2 = out.split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count -
                                                                      inter_tx_counter)].split('~')[1].split(':')[0]
            coord2_1 = out.split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count -
                                                                        inter_tx_counter)].split('~')[1].split(':')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom1) + ':' + str(coord1_1) + '~' + str(chrom2) + ':' + str(coord2_1)
            uniqname = read_name + '~' + bp_uname
            final.append(
                str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '~' +
                str(strands) + '\t' + str(chrom1) + '\t' + str('-') + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' +
                str(uniqname) + '\t' + str(signature)
            )
            inter_tx_counter = inter_tx_counter - 1
        else:
            raise Exception("Error: Unidentifiable SV")
    return final


# Generate SV signature
def getsignature(bp1, name, querymap, sign, realnmaps, elist, bitlist, bps_no, complex_sv, nchrom, short, piden, mismatch,
                 gap_ratio, strandness, sig_index):
    signlist = sign.split(',')
    if name == 'S-Nov_Ins':
        signature = ','.join(signlist[0:3]) + ',0,0,' + str(elist[0]) + ',0,' + str(bitlist[0]) + ',0,' + \
                    str(complex_sv) + ',' + str(realnmaps) + ',' + str(bps_no) + ',' + str(nchrom) + ',' + str(short) + ',' + \
                    str(piden[0]) + ',0,' + str(mismatch[0]) + ',0,' + str(gap_ratio[0]) + ',0'
        strands = str(strandness[0])
    elif name == 'E-Nov_Ins':
        signature = ','.join(signlist[-3:]) + ',0,0,' + str(elist[-1]) + ',0,' + str(bitlist[-1]) + ',0,' + str(complex_sv) + \
                    ',' + str(realnmaps) + ',' + str(bps_no) + ',' + str(nchrom) + ',' + str(short) + ',' + str(piden[-1]) + \
                    ',0,' + str(mismatch[-1]) + ',0,' + str(gap_ratio[-1]) + ',0'
        strands = str(strandness[-1])
    else:
        maps = ast.literal_eval(querymap)
        nmaps = len(maps)
        signature = ''
        strands = ''
        for i in range(nmaps):
            if int(maps[i][1]) == bp1:
                signature = ','.join(signlist[sig_index[i]:sig_index[i] + 5]) + ',' + str(elist[i]) + ',' + str(elist[i + 1]) + \
                            ',' + str(bitlist[i]) + ',' + str(bitlist[i + 1]) + ',' + str(complex_sv) + ',' + str(realnmaps) + \
                            ',' + str(bps_no) + ',' + str(nchrom) + ',' + str(short) + ',' + str(piden[i]) + ',' + \
                            str(piden[i + 1]) + ',' + str(mismatch[i]) + ',' + str(mismatch[i + 1]) + ',' + str(gap_ratio[i]) + \
                            ',' + str(gap_ratio[i + 1])
                strands = str(strandness[i]) + ',' + str(strandness[i + 1])
                break
    return signature, strands

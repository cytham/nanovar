"""
Parses BAM file and analyse CIGARS

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

import pysam
import os
import logging
from nv_detect_algo import sv_detect
from nv_parser import entry_parser, breakpoint_parser


def bam_parse(bam, unsigned int minlen, float splitpct, unsigned int minalign, str wk_dir, filter_file):
    cdef:
        unsigned int readlen, rstart, rend, flag, qlen, nm, seed
        unsigned long long basecov
        float ovlt
        int total_score
        str qname
        str rname
        list qseg, sseg, del_list, ins_list, cigar_tup, total_subdata, total_lines, contig_collect, total_out, sig_index
        bint adv
        object seg
        object sam = pysam.AlignmentFile(bam, "rb")
        dict repeat_dict = {}
        dict rlendict = {}
        dict main_dict = {}
        dict parse_dict = {}
        dict gapdict = {}
    seed = 0
    basecov = 0
    ovlt = 0.9  # Set overlap tolerance
    sig_index = [0, 2, 4]
    fasta = open(os.path.join(wk_dir, 'temp1.fa'), 'w')
    fasta2 = open(os.path.join(wk_dir, 'temp2.fa'), 'w')
    for seg in sam:
        flag = seg.flag
        qname = seg.query_name
        if flag == 4:
            fasta2.write('>' + qname + '\n' + seg.query_sequence + '\n')
            rlendict[qname] = len(seg.query_sequence)
            continue
        rname = seg.reference_name
        readlen = seg.infer_read_length()
        rstart = seg.reference_start
        rend = seg.reference_end
        qlen = seg.query_alignment_length
        nm = seg.get_tag('NM')
        total_score = seg.get_tag('AS')
        cigar_tup = seg.cigartuples
        adv, qseg, sseg, del_list, ins_list = read_cigar(cigar_tup, minlen, splitpct, rstart, rend, readlen)
        if flag in (0, 16):
            try:
                if repeat_dict[qname]:
                    pass
            except KeyError:
                repeat_dict[qname] = ''
                fasta.write(qname + '\n' + seg.query_sequence + '\n')
                rlendict[qname] = readlen
        try:
            main_dict[qname].append((adv, qname, rname, rstart, rend, readlen, qlen, flag, nm, total_score, qseg, sseg, del_list,
            ins_list))
        except KeyError:
            main_dict[qname] = []
            main_dict[qname].append((adv, qname, rname, rstart, rend, readlen, qlen, flag, nm, total_score, qseg, sseg, del_list,
            ins_list))
    fasta.close()
    fasta2.close()
    total_subdata, total_lines, contig_collect, total_out = [], [], [], []
    # Make gap dictionary
    gapdict = makegapdict(filter_file)
    for qname in main_dict:
        if len(main_dict[qname]) == 1:  # Single alignment read
            if not main_dict[qname][0][0]:  # if no sub-segments and not clipped read
                aln = main_dict[qname][0]
                qname = aln[1]
                rname = aln[2]
                rstart = aln[3]
                rend = aln[4]
                total_subdata.append(info_parse_simple(qname, rname, rstart, rend, minalign))
                # Add to base coverage
                basecov += rend - rstart
            else:  # if multiple sub-segments or clipped read
                seed += 1
                aln = main_dict[qname][0]
                qname = aln[1]
                rname = aln[2]
                readlen = aln[5]
                qlen = aln[6]
                flag = aln[7]
                nm = aln[8]
                total_score = aln[9]
                qseg = aln[10]
                sseg = aln[11]
                del_list = aln[12]
                ins_list = aln[13]
                total_lines, contig_collect = info_parse(qname, rname, readlen, qlen, flag, nm, total_score, qseg, sseg,
                del_list, ins_list, minalign)
                lines_sort = sorted(sorted(total_lines, key=lambda x: x[1], reverse=True), key=lambda y: y[0])
                temp1 = [tup[2] for tup in lines_sort]
                # Parse entries and correct overlap alignments
                subdata = entry_parser(temp1, contig_collect, ovlt)
                for entry in subdata:
                    total_subdata.append('\t'.join(entry.split('\t')[0:5]))
                    # Add to base coverage
                    basecov += int(entry.split('\t')[2])
                # SV detection
                out1, out2 = sv_detect(subdata, splitpct, minalign, gapdict)
                if out2 == '':
                    pass
                else:
                    # Parse breakpoints
                    final = breakpoint_parser(out2, minlen, sig_index, seed)
                    total_out.extend(final)
                    for i in final:
                        parse_dict[i.split('\t')[8]] = i
        else:  # Multiple alignment read
            seed += 1
            total_lines = []
            contig_collect = []
            for aln in main_dict[qname]:
                qname = aln[1]
                rname = aln[2]
                readlen = aln[5]
                qlen = aln[6]
                flag = aln[7]
                nm = aln[8]
                total_score = aln[9]
                qseg = aln[10]
                sseg = aln[11]
                del_list = aln[12]
                ins_list = aln[13]
                lines, contigs = info_parse(qname, rname, readlen, qlen, flag, nm, total_score, qseg, sseg,
                del_list, ins_list, minalign)
                total_lines.extend(lines)
                contig_collect.extend(contigs)
            contig_collect = sorted(set(contig_collect))
            lines_sort = sorted(sorted(total_lines, key=lambda x: x[1], reverse=True), key=lambda y: y[0])
            temp1 = [tup[2] for tup in lines_sort]
            # Parse entries and correct overlap alignments
            subdata = entry_parser(temp1, contig_collect, ovlt)
            for entry in subdata:
                total_subdata.append('\t'.join(entry.split('\t')[0:5]))
                # Add to base coverage
                basecov += int(entry.split('\t')[2])
            # SV detection
            out1, out2 = sv_detect(subdata, splitpct, minalign, gapdict)
            if out2 == '':
                pass
            else:
                # Parse breakpoints
                final = breakpoint_parser(out2, minlen, sig_index, seed)
                total_out.extend(final)
                for i in final:
                    parse_dict[i.split('\t')[8]] = i
    return total_subdata, total_out, basecov, parse_dict, rlendict, len(repeat_dict)


# Analyze CIGAR for Indels in segment and return read advancement call
cdef read_cigar(list cigar_tup, float minlen, float splitpct, unsigned int rstart, unsigned int rend, unsigned int readlen):
    cdef unsigned int qstart = 1
    cdef unsigned int softstart, softend, qlast, sstart, slast, dels, ins, b, v, t, move
    cdef list qseg = []
    cdef list sseg = []
    cdef list del_list = []
    cdef list ins_list = []

    if cigar_tup[0][0] == 4 or cigar_tup[0][0] == 5:
        qstart += cigar_tup[0][1]
        softstart = cigar_tup[0][1]
    else:
        softstart = 0
    if cigar_tup[-1][0] == 4 or cigar_tup[-1][0] == 5:
        softend = cigar_tup[-1][1]
    else:
        softend = 0
    qlast, dels, ins = 0, 0, 0
    sstart = rstart
    slast = sstart - 1
    for t, move in cigar_tup:
        if t in (4, 5):  # S or H
            qlast += move
        elif t == 0:  # M
            qlast += move
            slast += move
        elif t == 2:  # D
            dels += move
            if move > minlen:
                qseg.append([qstart, qlast])
                qstart = qlast + 1
                sseg.append([sstart, slast])
                slast += move
                sstart = slast + 1
                del_list.append(dels)
                ins_list.append(ins)
                dels = 0
                ins = 0
            else:
                slast += move
        elif t == 1:  # I
            ins += move
            if move > minlen:
                qseg.append([qstart, qlast])
                qlast += move
                qstart = qlast + 1
                sseg.append([sstart, slast])
                sstart = slast + 1
                del_list.append(dels)
                ins_list.append(ins)
                dels = 0
                ins = 0
            else:
                qlast += move
        else:
            logging.critical('Error: Unrecognized CIGAR translated symbol "%s"' % t)
            raise Exception('Error: Unrecognized CIGAR translated symbol "%s"' % t)
    # Last segment
    qseg.append([qstart, qlast - softend])
    sseg.append([sstart, slast])
    del_list.append(dels)
    ins_list.append(ins)
    if len(qseg) == 1:
        if softstart > splitpct * readlen or softend > splitpct * readlen:
            return True, qseg, sseg, del_list, ins_list
        else:
            return False, qseg, sseg, del_list, ins_list
    else:
        return True, qseg, sseg, del_list, ins_list



def info_parse(qname, rname, readlen, qlen, flag, nm, total_score, qseg, sseg, del_list, ins_list, minlen=200):
    lines = []
    total_del = sum([x for x in del_list])
    total_ins = sum([x for x in ins_list])
    total_gap = total_del + total_ins
    total_mismatch = nm - total_gap
    nsegs = len(qseg)
    contig_collect = []
    for i in range(nsegs):
        substart = sseg[i][0]
        substretch = sseg[i][1] - sseg[i][0]
        if substretch < minlen:
            continue
        qstart, qstretch, strand = query_sign(qseg[i][0], qseg[i][1], flag, readlen)
        ndel = del_list[i]
        nins = ins_list[i]
        gaps = ndel + nins
        subratio = round(qstretch / qlen, 3)
        mismatch = int(round(total_mismatch * subratio, 0))
        matches = max(qstretch - mismatch - nins, 1)
        pident = min(round(matches * 100 / substretch, 2), 100.00)
        score = int(round(total_score * subratio, 0))
        priority = align_priority(flag)
        line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rname, substart, substretch, '+',
                                                                           qname, qstart, qstretch, strand,
                                                                           readlen, 0.0, score, pident, mismatch, gaps)
        contig_collect.append(rname)
        lines.append((priority, score, line))
    return lines, contig_collect

# Parse segments without breakends
def info_parse_simple(qname, rname, substart, rend, minlen=200):
    substretch = rend - substart
    if substretch < minlen:
        pass
    line = "%s\t%s\t%s\t%s\t%s" % (rname, substart, substretch, '+', qname)
    return line

# Parse query start and end values according to read direction
def query_sign(start, end, flag, readlen):
    if flag in (0, 2048, 256):
        return start, end - start, '+'
    elif flag in (16, 2064, 272):
        return readlen - end, end - start, '-'
    else:
        return start, end - start, '+'

# Set alignment priority according to flag
def align_priority(flag):
    if flag in (0, 16):
        return 1
    elif flag in (2048, 2064):
        return 2
    else:
        return 3


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

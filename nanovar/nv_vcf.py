"""
Functions for creating VCF file.

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
import math
from datetime import date
from natsort import natsorted
from nanovar import __version__


def create_vcf(wk_dir, thres, nn_out, ref_path, read_path, read_name, mm_cmd, contig_len_dict, homo_t, het_t, minlen, depth,
               index2te, nv_cmd, alt_seq):
    rdata = nn_out
    # Calculating number of entries
    t = len(rdata)
    # Set range according to number of queries
    k = range(t)
    # Create last line dummy
    rdata.append('dum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum')
    ref_path = os.path.abspath(ref_path)
    read_path = os.path.abspath(read_path)
    vcf_total = open(os.path.join(wk_dir, '%s.nanovar.total.vcf' % read_name), 'w')
    vcf_pass = open(os.path.join(wk_dir, '%s.nanovar.pass.vcf' % read_name), 'w')
    add_header(vcf_total, read_path, ref_path, mm_cmd, read_name, contig_len_dict, thres, depth, nv_cmd)
    add_header(vcf_pass, read_path, ref_path, mm_cmd, read_name, contig_len_dict, thres, depth, nv_cmd)
    sv_sup = open(os.path.join(wk_dir, 'sv_support_reads.tsv'), 'w')
    sv_sup.write('SV-ID\tSupporting_reads (readname~index1,readname~index2...)\n')
    tmpread = []
    out = []
    for i in k:
        tmpread.append(rdata[i])
        covl = int(tmpread[0].split('\t')[10])
        normcov = int(float(tmpread[0].split('\t')[12]))
        dnnscore = float(tmpread[0].split('\t')[13])
        dnn = str(round(float(dnnscore), 3))
        phred = phredc(dnnscore)
        sv_id = tmpread[0].split('\t')[6].split('~')[0]
        bp_name = tmpread[0].split('\t')[3].split(' ')[0]
        chrm1 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
        filt = filterer(float(phred), thres)
        phred = str(abs(round(phred, 1)))
        dp = str(covl + normcov)
        ratio = round(float(covl) / int(dp), 3)
        geno = genotyper(ratio, homo_t, het_t)
        if bp_name == 'Nov_Ins':
            sv = '<INS>'
            new_sv_id = 'NV.INS.' + sv_id.split('nv_')[1]
            sv_len = tmpread[0].split('\t')[3].split(' ')[1].split('~')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            if sv_id in index2te:
                te = ";TE=" + ','.join(sorted(index2te[sv_id]))
            else:
                te = ''
            ref_bases = alt_seq[sv_id][0]
            alt_bases = alt_seq[sv_id]
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\t{}\t{}\t'.format(ref_bases, alt_bases) + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=INS;END=' + str(coord1 + 1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(
                covl) + ';NN=' +
                       str(dnn) + te + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'E-Nov_Ins_bp' or bp_name == 'S-Nov_Ins_bp':
            pass
            sv = '<INS>'
            new_sv_id = 'NV.INS.' + sv_id.split('nv_')[1]
            sv_len = tmpread[0].split('\t')[3].split(' ')[1].split('~')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            if sv_id in index2te:
                te = ";TE=" + ','.join(sorted(index2te[sv_id]))
            else:
                te = ''
            ref_bases = alt_seq[sv_id][0]
            alt_bases = alt_seq[sv_id]
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\t{}\t{}\t'.format(ref_bases, alt_bases) + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=INS;END=' + str(coord1 + 1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(
                covl) + ';NN=' +
                       str(dnn) + te + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'Del':
            sv = '<DEL>'
            new_sv_id = 'NV.DEL.' + sv_id.split('nv_')[1]
            sv_len = int(tmpread[0].split('\t')[3].split(' ')[1].split('~')[0])
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            if coord2 - coord1 < minlen:
                mid = (coord2 + coord1) / 2
                coord1 = max(1, int(mid - round(sv_len/2, 0)))
                coord2 = int(mid + round(sv_len/2, 0) + 1)
            sv_len = '-' + str(coord2 - coord1)
            ref_bases = alt_seq[sv_id]
            alt_bases = alt_seq[sv_id][0]
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\t{}\t{}\t'.format(ref_bases, alt_bases) + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=DEL;END=' + str(coord2) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'Inv':
            sv = '<INV>'
            new_sv_id = 'NV.INV.' + sv_id.split('nv_')[1]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            if coord2 - coord1 < minlen:
                mid = (coord2 + coord1) / 2
                coord1 = max(1, int(mid - round(minlen / 2, 0)))
                coord2 = int(mid + round(minlen / 2, 0) + 1)
            sv_len = str(coord2 - coord1)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(sv) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=INV;END=' + str(coord2) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'Inv(1)' or bp_name == 'Inv(2)':
            sv = '<INV>'
            new_sv_id = 'NV.INV.' + sv_id.split('nv_')[1]
            # sv_len = tmpread[0].split('\t')[3].split(' ')[1].split('~')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            if coord2 - coord1 < minlen:
                mid = (coord2 + coord1) / 2
                coord1 = max(1, int(mid - round(minlen / 2, 0)))
                coord2 = int(mid + round(minlen / 2, 0) + 1)
            sv_len = str(coord2 - coord1)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(sv) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=INV;END=' + str(coord2) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'TDupl':
            sv = '<DUP>'
            new_sv_id = 'NV.DUP.' + sv_id.split('nv_')[1]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            if coord2 - coord1 < minlen:
                mid = (coord2 + coord1) / 2
                coord1 = max(1, int(mid - round(minlen / 2, 0)))
                coord2 = int(mid + round(minlen / 2, 0) + 1)
            sv_len = str(coord2 - coord1)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(sv) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=DUP;END=' + str(coord2) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'Intra-Ins':
            new_sv_id = 'NV.BND.' + sv_id.split('nv_')[1]
            sv_len = '0'
            strands = tmpread[0].split('\t')[3].split(' ')[1].split('~')[1].split(',')
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            a, b = breakend_alt(strands, chrm1, coord1, chrm2, coord2)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + ';SV2=TRA' + ';CHR2=' + str(chrm2) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TRA' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm2) + '\t' + str(coord2) + '\t' + str(sv_id) + '\tN\t' + str(b) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord2+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TRA' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'Intra-Ins(1)' or bp_name == 'Intra-Ins(2)':
            new_sv_id = 'NV.BND.' + sv_id.split('nv_')[1]
            sv_len = tmpread[0].split('\t')[3].split(' ')[1].split('~')[0]
            strands = tmpread[0].split('\t')[3].split(' ')[1].split('~')[1].split(',')
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            a, b = breakend_alt(strands, chrm1, coord1, chrm2, coord2)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + ';SV2=TPO' + ';CHR2=' + str(chrm2) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TPO' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm2) + '\t' + str(coord2) + '\t' + str(sv_id) + '\tN\t' + str(b) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord2+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TPO' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'Inter-Ins(1)' or bp_name == 'Inter-Ins(2)':
            new_sv_id = 'NV.BND.' + sv_id.split('nv_')[1]
            sv_len = tmpread[0].split('\t')[3].split(' ')[1].split('~')[0]
            strands = tmpread[0].split('\t')[3].split(' ')[1].split('~')[1].split(',')
            chrm2 = tmpread[0].split('\t')[6].split('~')[2].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[2].split(':')[1])
            a, b = breakend_alt(strands, chrm1, coord1, chrm2, coord2)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + ';SV2=TPO' + ';CHR2=' + str(chrm2) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TPO' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm2) + '\t' + str(coord2) + '\t' + str(sv_id) + '\tN\t' + str(b) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord2+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TPO' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        elif bp_name == 'InterTx':
            new_sv_id = 'NV.BND.' + sv_id.split('nv_')[1]
            sv_len = '0'
            strands = tmpread[0].split('\t')[3].split(' ')[1].split('~')[1].split(',')
            chrm2 = tmpread[0].split('\t')[6].split('~')[2].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[2].split(':')[1])
            a, b = breakend_alt(strands, chrm1, coord1, chrm2, coord2)
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(new_sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
                       filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
                       str(dnn) + ';SV2=TRA' + ';CHR2=' + str(chrm2) + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\tN\t' + str(a) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord1+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TRA' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
            # out.append(str(chrm2) + '\t' + str(coord2) + '\t' + str(sv_id) + '\tN\t' + str(b) + '\t' + str(phred) + '\t' +
            #            filt + '\t' + 'SVTYPE=BND;END=' + str(coord2+1) + ';SVLEN=' + str(sv_len) + ';SR=' + str(covl) + ';NN=' +
            #            str(dnn) + ';SV2=TRA' + '\tGT:DP:AD\t' + geno + ':' + dp + ':' + str(normcov) + ',' + str(covl))
        else:
            raise Exception("Error: Unrecognised breakpoint name")
        sv_sup.write(new_sv_id + '\t' + tmpread[0].split('\t')[8] + ',' + tmpread[0].split('\t')[11] + '\n')
        tmpread = []
    out_sort1 = natsorted(out, key=lambda x: x.split('\t')[1])  # sort by POS
    total_vcf = natsorted(out_sort1, key=lambda x: x.split('\t')[0])  # sort by CHROM
    for e in total_vcf:
        vcf_total.write(e + '\n')
        if e.split('\t')[6] == 'PASS':
            vcf_pass.write(e + '\n')
    vcf_total.close()
    vcf_pass.close()
    sv_sup.close()


# Write VCF header
def add_header(vcf_file, read_path, ref_path, mm_cmd, read_name, contig_len_dict, thres, depth, nv_cmd):
    today = date.today()
    today_date = today.strftime("%d-%m-%Y")
    vcf_file.write('##fileformat=VCFv4.2\n')
    vcf_file.write('##fileDate=%s\n' % today_date)
    vcf_file.write('##source=NanoVar-%s\n' % __version__)
    vcf_file.write('##source_reads=%s\n' % read_path)
    vcf_file.write('##reference=%s\n' % ref_path)
    vcf_file.write('##command=%s\n' % nv_cmd)
    vcf_file.write('##mapping=%s\n' % mm_cmd)
    vcf_file.write('##phasing=none\n')
    vcf_file.write('##depth_of_coverage=%sX\n' % str(depth))
    for key in contig_len_dict:
        vcf_file.write('##contig=<ID=%s,length=%s>\n' % (key, str(contig_len_dict[key])))
    vcf_file.write('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">\n')
    vcf_file.write('##ALT=<ID=DEL,Description="Deletion relative to the reference">\n')
    vcf_file.write('##ALT=<ID=INV,Description="Inversion relative to the reference">\n')
    vcf_file.write('##ALT=<ID=DUP,Description="Tandem duplication relative to the reference">\n')
    vcf_file.write('##ALT=<ID=BND,Description="Breakend relative to the reference">\n')
    vcf_file.write('##FILTER=<ID=PASS,Description="Quality above/equal to %s">\n' % str(thres))
    vcf_file.write('##FILTER=<ID=FAIL,Description="Quality below %s">\n' % str(thres))
    vcf_file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    vcf_file.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the described variant">\n')
    vcf_file.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Estimated length of the variant">\n')
    vcf_file.write('##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of supporting reads">\n')
    vcf_file.write('##INFO=<ID=NN,Number=1,Type=Float,Description="Neural network confidence probability">\n')
    vcf_file.write('##INFO=<ID=SV2,Number=1,Type=String,Description="BND SV assessment: TPO - Transposition, '
                   'TRA - Translocation">\n')
    vcf_file.write('##INFO=<ID=TE,Number=.,Type=String,Description="Identity of transposable element insertion">\n')
    vcf_file.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Second chromosome name for BND SVs">\n')
    vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    vcf_file.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">\n')
    vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % read_name)


# Phred quality calculator
def phredc(p):
    probfalse = 1.000000001 - float(p)
    return math.log10(probfalse)*(-10)


# Filter based on threshold score
def filterer(x, thres):
    if x < thres:
        return "FAIL"
    else:
        return "PASS"


# Predict SV genotype
def genotyper(x, homo_t, het_t):
    if x >= homo_t:
        return "1/1"
    elif x >= het_t:
        return "0/1"
    else:
        return "./."


# Construct breakend alteration
def breakend_alt(strand_list, chrm1, coord1, chrm2, coord2):
    coord1 = str(coord1)
    coord2 = str(coord2)
    if strand_list == ['+', '+']:
        a = 'N[' + chrm2 + ':' + coord2 + '['
        b = ']' + chrm1 + ':' + coord1 + ']N'
    elif strand_list == ['+', '-']:
        a = 'N]' + chrm2 + ':' + coord2 + ']'
        b = 'N]' + chrm1 + ':' + coord1 + ']'
    elif strand_list == ['-', '+']:
        a = '[' + chrm2 + ':' + coord2 + '[N'
        b = '[' + chrm1 + ':' + coord1 + '[N'
    elif strand_list == ['-', '-']:
        a = ']' + chrm2 + ':' + coord2 + ']N'
        b = 'N[' + chrm1 + ':' + coord1 + '['
    else:
        raise Exception('Strand information not found.')
    return a, b

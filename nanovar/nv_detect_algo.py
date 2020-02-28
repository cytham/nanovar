"""
SV detection algorithm.

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


def sv_detect(subdata, splitpct, minalign, gapdict):
    c = len(subdata)  # Number of entries
    if c == 0:
        out1, out2 = '', ''
    else:
        if c == 1:
            sv_status = single_alignment_check(subdata, splitpct, minalign)
        else:
            subdata = reduce_multi_map(subdata)
            c = len(subdata)  # Recalculate number of entries
            sv_status = True
        if sv_status:
            tmpread = subdata
            actual_no_maps = subdata[0].split('\t')[14].split('=')[1]
            chromocollect = tmpread[0].split('\t')[15].split('c')[0]
            # Set minimum length for read terminal breakends
            termlength = 100
            # Temporary set SV len cutoff
            ncutoff = 19
            j = range(c)
            p = range(c - 1)
            # Define SV flags
            strand_total, evalue_total, bitscore_ratio_total, piden_total, mismatch_ratio_total, gap_ratio_total, query, \
                subject, chromorder, qurygaps, complex_sv, iso_nov_ins, iso_ins_size, iso_ins_range, sv_range, intra_ins, \
                intra_ins_range, dele, del_size, del_range, tdup, tdup_range, nov_ins, ins_size, ins_range, inv, inv_range, \
                complex_nov_ins, complex_ins_size, inter_ins, inter_ins_range, intertx, intertx_break = ([] for i in range(33))
            # Satellite elements counter, Short alignments counter, Ref chromosome counter, Strandness counter, Second bp markers
            a = f = q = s = g = gtx = 0
            # Choose tentative reference chromosome and strandness
            refchr = tmpread[0].split('\t')[0]
            refstrand = tmpread[0].split('\t')[7]
            # Retrieving read length
            readlength = int(tmpread[0].split('\t')[8])
            # Define short length alignment (Arbitrary definition as <5% of read length)
            t = int(readlength) * 0.05
            for i in j:  # Collecting alignment information and filtering
                # Scale evalue 100000*, round to 5 dec places
                evalue = round(float(tmpread[i].split('\t')[9]) * 100000, 5)
                # Calculate bitscore ratio
                bitscore_ratio = round(float(tmpread[i].split('\t')[10]) / int(tmpread[i].split('\t')[6]), 5)
                # Collect alignment percentage identity
                piden = tmpread[i].split('\t')[11]
                # Calculate mismatch ratio
                mismatch_ratio = round(float(tmpread[i].split('\t')[12]) / int(tmpread[i].split('\t')[6]), 5)
                # Calculate gap ratio
                gap_ratio = round(float(tmpread[i].split('\t')[13]) / int(tmpread[i].split('\t')[6]), 5)
                # Get strandness
                strand = tmpread[i].split('\t')[7]
                # Collect all alignment info
                evalue_total.append(str(evalue))
                bitscore_ratio_total.append(str(bitscore_ratio))
                piden_total.append(str(piden))
                mismatch_ratio_total.append(str(mismatch_ratio))
                gap_ratio_total.append(str(gap_ratio))
                strand_total.append(str(strand))
                # Collecting query and subject information
                query.append([int(tmpread[i].split('\t')[5]), int(tmpread[i].split('\t')[5]) + int(tmpread[i].split('\t')[6])])
                if tmpread[i].split('\t')[7] == '+':
                    subject.append([int(tmpread[i].split('\t')[1]), int(tmpread[i].split('\t')[1]) +
                                    int(tmpread[i].split('\t')[2])])
                elif tmpread[i].split('\t')[7] == '-':
                    subject.append([int(tmpread[i].split('\t')[1]) + int(tmpread[i].split('\t')[2]),
                                    int(tmpread[i].split('\t')[1])])
                # Filtering out strict telomeric and centromere regions
                chrm = tmpread[i].split('\t')[0].strip()
                chromorder.append(chrm)
                subjstart = int(tmpread[i].split('\t')[1])
                subjend = int(tmpread[i].split('\t')[1]) + int(tmpread[i].split('\t')[2])
                if gapdict is not None:
                    try:
                        for v in gapdict[chrm]:
                            if v[0] <= (subjstart or subjend) <= v[1]:
                                a += 1
                    except KeyError:
                        pass
                # Count short alignments (Arbitrary definition as <5% of read length)
                if int(tmpread[i].split('\t')[6]) < t:
                    f += 1
                # Count support for ref chromosome and strandness
                if tmpread[i].split('\t')[0] == refchr:
                    q += 1
                if tmpread[i].split('\t')[7] == refstrand:
                    s += 1
            # Calculating query middle gaps length
            for i in p:
                qurygaps.append(query[i + 1][0] - query[i][1])
            # Check for Satellite element
            if a == c:  # All alignments lies within satellite region
                sv_status = False
            elif 0 < a:
                sv_status = False
            # Check for short alignment lengths
            short = float(f) / c
            # Defining major strandness
            majstrand = ''  # Dummy
            if float(s) / c >= 0.5:
                majstrand = refstrand
            else:
                if refstrand == '+':
                    majstrand = '-'
                elif refstrand == '-':
                    majstrand = '+'
            # Counting number of inversions
            if majstrand == refstrand:
                if (c - s) > 1:
                    complex_sv.append('Multi_inv')
            else:
                if s > 1:
                    complex_sv.append('Multi_inv')
            # Calculating start/end gaps length and evaluation based on 100bp or more than ncutoff bases
            startgap = int(tmpread[0].split('\t')[5])
            endgap = int(tmpread[int(c) - 1].split('\t')[8]) - int(tmpread[int(c) - 1].split('\t')[6]) - int(
                tmpread[int(c) - 1].split('\t')[5])
            if startgap > termlength and startgap > ncutoff:
                iso_nov_ins.append('S-Nov-Ins')
                iso_ins_size.append('S-Nov-Ins' + str(startgap))
                iso_ins_range.append(str(tmpread[0].split('\t')[0].strip()) + ':' + str(subject[0][0]))
                sv_range.append('S-Nov_Ins:1-' + str(query[0][0] + 1))
            if endgap > termlength and endgap > ncutoff:
                iso_nov_ins.append('E-Nov-Ins')
                iso_ins_size.append('E-Nov-Ins' + str(endgap))
                iso_ins_range.append(str(tmpread[-1].split('\t')[0].strip()) + ':' + str(subject[-1][1]))
                sv_range.append('E-Nov_Ins:' + str(query[-1][1] - 1) + '-' + str(readlength))
            # Begin scanning all alignments on read
            if sv_status:
                for i in p:
                    sbp = 0
                    qurygap = query[i + 1][0] - query[i][1]
                    # Same chromosome
                    if tmpread[i].split('\t')[0] == tmpread[i + 1].split('\t')[0]:
                        # Same strandness (1)
                        if tmpread[i].split('\t')[7] == '+' and tmpread[i + 1].split('\t')[7] == '+':
                            if i != int(g - 1):  # Prevent duplicate record
                                # Gapstudy
                                subjgap = subject[i + 1][0] - subject[i][1]
                                if subjgap > ncutoff:  # If gap is positive and more than ncutoff
                                    # Subject gap needs to be at least twice the size of query gap
                                    if float(subjgap) / (qurygap + 1) >= 2:
                                        # First breakpoint detected
                                        # Detecting second breakpoint
                                        for u in range(i + 2, c):
                                            size = query[u][0] - query[i][1]
                                            if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                                    tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                                    any(subject[i][1] < y < min(subject[i+1][0], subject[i][1]+size+1000) for y
                                                        in subject[u]):
                                                # Second breakpoint detected, insertion event
                                                sbp = 1
                                                g = int(u)
                                                break
                                        if sbp == 1:
                                            # Second breakpoint present
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                            # Duplicate record for Intra-Ins(1)
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                                   str(subject[i][1]) + '-' + str(subject[i + 1][0]))
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                                   str(subject[g - 1][1]) + '-' + str(subject[g][0]))
                                            sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) +
                                                            ',' + 'Intra-Ins(2):' + str(query[g - 1][1]) + '-' +
                                                            str(query[g][0]))
                                        elif sbp == 0:
                                            dele.append('Del')
                                            del_size.append('Del' + str(subjgap-qurygap))
                                            del_range.append(
                                                str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(
                                                    subject[i + 1][0]))
                                            sv_range.append('Del:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If gap is negative and smaller than 100000 (Arbitrary tandemDup size gauge)
                                elif -100000 < subjgap < -20:
                                    tdup.append('TDupl')
                                    tdup_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                      str(subject[i + 1][0]))
                                    sv_range.append('TDupl:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If gap is negative and larger than or equal to 100000
                                elif subjgap <= -100000:
                                    # Detecting second breakpoint
                                    for u in range(i + 2, c):
                                        size = query[u][0] - query[i][1]
                                        if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                                tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                                any(subject[i][1] < y < subject[i][1]+size+1000 for y in subject[u]):  # May
                                            # still contain large deletion
                                            sbp = 1
                                            g = int(u)
                                            break
                                    if sbp == 1:
                                        # Second breakpoint present
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                        # Duplicate record for Intra-Ins(1)
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                               str(subject[i][1]) + '-' + str(subject[i + 1][0]))
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                               str(subject[g - 1][1]) + '-' + str(subject[g][0]))
                                        sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) + ',' +
                                                        'Intra-Ins(2):' + str(query[g - 1][1]) + '-' + str(query[g][0]))
                                    elif sbp == 0:
                                        intra_ins.append('Intra-Ins')
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                               str(subject[i][1]) + '-' + str(subject[i + 1][0]))
                                        sv_range.append('Intra-Ins:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If query gap more than cutoff and subject gap is positive
                                if qurygap > ncutoff and subjgap >= 0:
                                    if float(qurygap) / (subjgap + 1) >= 2:  # Test for novel insertion ratio
                                        nov_ins.append('Nov-Ins')
                                        ins_size.append('Nov-Ins' + str(qurygap-subjgap))
                                        ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) +
                                                         '-' + str(int(subject[i][1]) + 1))
                                        sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If query gap more than cutoff and subject gap is negative
                                elif qurygap > ncutoff and subjgap < 0:
                                    nov_ins.append('Nov-Ins')
                                    ins_size.append('Nov-Ins' + str(qurygap-subjgap))
                                    ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                     str(int(subject[i][1]) + 1))
                                    sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                        # Same strandness (2)
                        if tmpread[i].split('\t')[7] == '-' and tmpread[i + 1].split('\t')[7] == '-':
                            if i != int(g - 1):  # Prevent duplicate record
                                # Gapstudy
                                subjgap = subject[i][1] - subject[i + 1][0]
                                if subjgap > ncutoff:  # Is gap positive value and more than ncutoff?
                                    if float(subjgap) / (qurygap + 1) >= 2:  # Test for significant gap using ratio
                                        # First breakpoint detected
                                        # Detecting second breakpoint
                                        for u in range(i + 2, c):
                                            size = query[u][0] - query[i][1]
                                            if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                                    tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                                    any(max(subject[i + 1][0], subject[i][1]-size-1000) < y < subject[i][1] for
                                                        y in subject[u]):
                                                # Second breakpoint detected, insertion event
                                                sbp = 1
                                                g = int(u)
                                                break
                                        if sbp == 1:
                                            # Second breakpoint present
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                            # Duplicate record for Intra-Ins(1)
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                                   str(subject[i][1]) + '-' + str(subject[i + 1][0]))
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                                   str(subject[g - 1][1]) + '-' + str(subject[g][0]))
                                            sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) +
                                                            ',' + 'Intra-Ins(2):' + str(query[g - 1][1]) + '-' + str(query[g][0]))
                                        elif sbp == 0:
                                            dele.append('Del')
                                            del_size.append('Del' + str(subjgap-qurygap))
                                            del_range.append(
                                                str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(
                                                    subject[i + 1][0]))
                                            sv_range.append('Del:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If gap is negative and smaller than 100000 (Arbitrary tandemDup size gauge)
                                elif -100000 < subjgap < -20:
                                    tdup.append('TDupl')
                                    tdup_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                      str(subject[i + 1][0]))
                                    sv_range.append('TDupl:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If gap is negative and larger than or equal to 100000
                                elif subjgap <= -100000:
                                    # Detecting second breakpoint
                                    for u in range(i + 2, c):
                                        size = query[u][0] - query[i][1]
                                        if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                                tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                                any(subject[i][1]-size-1000 < y < subject[i][1] for y in subject[u]):
                                            # May still contain large deletion
                                            sbp = 1
                                            g = int(u)
                                            break
                                    if sbp == 1:
                                        # Second breakpoint present
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                        # Duplicate record for Intra-Ins(1)
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0] - query[i][1]))
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                               str(subject[i][1]) + '-' + str(subject[i + 1][0]))
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                               str(subject[g - 1][1]) + '-' + str(subject[g][0]))
                                        sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) + ',' +
                                                        'Intra-Ins(2):' + str(query[g - 1][1]) + '-' + str(query[g][0]))
                                    elif sbp == 0:
                                        intra_ins.append('Intra-Ins')
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' +
                                                               str(subject[i][1]) + '-' + str(subject[i + 1][0]))
                                        sv_range.append('Intra-Ins:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If query gap more than cutoff and subject gap is positive
                                if qurygap > ncutoff and subjgap >= 0:
                                    if float(qurygap) / (subjgap + 1) >= 2:  # Test for novel insertion ratio
                                        nov_ins.append('Nov-Ins')
                                        ins_size.append('Nov-Ins' + str(qurygap-subjgap))
                                        ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                         str(int(subject[i][1]) + 1))
                                        sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                # If query gap more than cutoff and subject gap is negative
                                elif qurygap > ncutoff and subjgap < 0:
                                    nov_ins.append('Nov-Ins')
                                    ins_size.append('Nov-Ins' + str(qurygap-subjgap))
                                    ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                     str(int(subject[i][1]) + 1))
                                    sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                        # Different strandness (1)
                        elif tmpread[i].split('\t')[7] == '+' and tmpread[i + 1].split('\t')[7] == '-':  # Inversion
                            if i != int(g - 1):
                                for u in range(i + 2, c):
                                    size = query[u][0] - query[i][1]
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                            tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                            any(subject[i][1] < y < subject[i][1]+size+1000 for y in subject[u]):
                                        # Second breakpoint detected, paired inversion event
                                        sbp = 1
                                        g = int(u)
                                        break
                                if sbp == 1:
                                    # Second breakpoint present
                                    inv.append('Pair_Inv' + str(query[g][0] - query[i][1]))
                                    # Duplicate record for Inv(1)
                                    inv.append('Pair_Inv' + str(query[g][0] - query[i][1]))
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                     str(subject[i + 1][0]))
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g - 1][1]) + '-' +
                                                     str(subject[g][0]))
                                    sv_range.append('Inv(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) + ',' + 'Inv(2):' +
                                                    str(query[g - 1][1]) + '-' + str(query[g][0]))
                                elif sbp == 0:
                                    inv.append('Unpair_Inv')
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                     str(subject[i + 1][0]))
                                    sv_range.append('Inv:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                            if qurygap > ncutoff:  # If query gap more than ncutoff
                                complex_nov_ins.append('Inv-Nov-Ins')
                                complex_sv.append('Inv-Nov_Ins')
                                complex_ins_size.append('Inv-Nov_Ins' + str(qurygap))
                        # Different Strandness(2)
                        elif tmpread[i].split('\t')[7] == '-' and tmpread[i + 1].split('\t')[7] == '+':  # Inversion
                            if i != int(g - 1):
                                for u in range(i + 2, c):
                                    size = query[u][0] - query[i][1]
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                            tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                            any(subject[i][1]-size-1000 < y < subject[i][1] for y in subject[u]):
                                        # Second breakpoint detected, paired inversion event
                                        sbp = 1
                                        g = int(u)
                                        break
                                if sbp == 1:
                                    # Second breakpoint present
                                    inv.append('Pair_Inv' + str(query[g][0] - query[i][1]))
                                    # Duplicate record for Inv(1)
                                    inv.append('Pair_Inv' + str(query[g][0] - query[i][1]))
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                     str(subject[i + 1][0]))
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g - 1][1]) + '-' +
                                                     str(subject[g][0]))
                                    sv_range.append('Inv(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) + ',' + 'Inv(2):' +
                                                    str(query[g - 1][1]) + '-' + str(query[g][0]))
                                elif sbp == 0:
                                    inv.append('Unpair_Inv')
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' +
                                                     str(subject[i + 1][0]))
                                    sv_range.append('Inv:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                            if qurygap > ncutoff:  # If query gap more than ncutoff
                                complex_nov_ins.append('Inv-Nov-Ins')
                                complex_sv.append('Inv-Nov_Ins')
                                complex_ins_size.append('Inv-Nov_Ins' + str(qurygap))
                    # Different chromosome
                    elif tmpread[i].split('\t')[0] != tmpread[i + 1].split('\t')[0]:
                        # Strandness+
                        if tmpread[i].split('\t')[7] == '+':
                            if i != int(gtx - 1):
                                for u in range(i + 2, c):
                                    size = query[u][0] - query[i][1]
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                            tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                            any(subject[i][1] < y < subject[i][1]+size+1000 for y in subject[u]):  # Detecting
                                        # second breakpoint
                                        # Second breakpoint detected, insertion event
                                        sbp = 1
                                        gtx = int(u)
                                        break
                                if sbp == 1:
                                    # Second breakpoint present
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0] - query[i][1]))
                                    # Duplicate record for Inter-Ins(1)
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0] - query[i][1]))
                                    inter_ins_range.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' +
                                                           str(tmpread[i + 1].split('\t')[0]) + ':' + str(subject[i + 1][0]))
                                    inter_ins_range.append(str(tmpread[gtx - 1].split('\t')[0]) + ':' +
                                                           str(subject[gtx - 1][1]) + '~' + str(tmpread[gtx].split('\t')[0]) +
                                                           ':' + str(subject[gtx][0]))
                                    sv_range.append('Inter-Ins(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) + ',' +
                                                    'Inter-Ins(2):' + str(query[gtx - 1][1]) + '-' + str(query[gtx][0]))
                                elif sbp == 0:
                                    intertx.append('InterTx')
                                    intertx_break.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' +
                                                         str(tmpread[i + 1].split('\t')[0]) + ':' + str(subject[i + 1][0]))
                                    sv_range.append('InterTx:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                if qurygap > ncutoff:  # If query gap more than ncutoff
                                    complex_nov_ins.append('Inter-Nov-Ins')
                                    complex_sv.append('Inter-Nov_Ins')
                                    complex_ins_size.append('Inter-Nov_Ins' + str(qurygap))
                        # Strandness-
                        if tmpread[i].split('\t')[7] == '-':
                            if i != int(gtx - 1):
                                for u in range(i + 2, c):
                                    size = query[u][0] - query[i][1]
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and \
                                            tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and \
                                            any(subject[i][1]-size-1000 < y < subject[i][1] for y in subject[u]):  # Detecting
                                        # second breakpoint
                                        # Second breakpoint detected, insertion event
                                        sbp = 1
                                        gtx = int(u)
                                        break
                                if sbp == 1:
                                    # Second breakpoint present
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0] - query[i][1]))
                                    # Duplicate record for Inter-Ins(1)
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0] - query[i][1]))
                                    inter_ins_range.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' +
                                                           str(tmpread[i + 1].split('\t')[0]) + ':' + str(subject[i + 1][0]))
                                    inter_ins_range.append(str(tmpread[gtx - 1].split('\t')[0]) + ':' +
                                                           str(subject[gtx - 1][1]) + '~' + str(tmpread[gtx].split('\t')[0]) +
                                                           ':' + str(subject[gtx][0]))
                                    sv_range.append('Inter-Ins(1):' + str(query[i][1]) + '-' + str(query[i + 1][0]) + ',' +
                                                    'Inter-Ins(2):' + str(query[gtx - 1][1]) + '-' + str(query[gtx][0]))
                                elif sbp == 0:
                                    intertx.append('InterTx')
                                    intertx_break.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' + str(
                                        tmpread[i + 1].split('\t')[0]) + ':' + str(subject[i + 1][0]))
                                    sv_range.append('InterTx:' + str(query[i][1]) + '-' + str(query[i + 1][0]))
                                if qurygap > ncutoff:  # If query gap more than ncutoff
                                    complex_nov_ins.append('Inter-Nov-Ins')
                                    complex_sv.append('Inter-Nov_Ins')
                                    complex_ins_size.append('Inter-Nov_Ins' + str(qurygap))
            if dele == [] and iso_nov_ins == [] and nov_ins == [] and complex_nov_ins == [] and tdup == [] and inv == [] and \
                    intra_ins == [] and inter_ins == [] and intertx == [] and complex_sv == []:
                sv_status = False
            total_sv_complex = ' '.join(dele) + ' ' + ' '.join(iso_nov_ins) + ' ' + ' '.join(nov_ins) + ' ' + ' '.join(
                tdup) + ' ' + ' '.join(inv) + ' ' + ' '.join(intra_ins) + ' ' + ' '.join(inter_ins) + ' ' + ' '.join(intertx)
            if sv_status:
                querypercent = []
                qurygappercent = []
                fullquerypercent = []
                for i in p:  # Calculating gap percentages
                    qurygappercent.append(str((float(qurygaps[i]) / readlength) * 100))
                startgappercent = str((float(startgap) / readlength) * 100)
                qurygappercent.append(str((float(endgap) / readlength) * 100))
                for i in j:
                    querypercent.append(
                        str((float(tmpread[i].split('\t')[6]) / readlength) * 100))  # Calculating query map percentage
                fullquerypercent.append('%.1f' % float(startgappercent) + '%')
                for i in j:
                    fullquerypercent.append('(' + '%.1f' % float(querypercent[i]) + '%)')
                    fullquerypercent.append('%.1f' % float(qurygappercent[i]) + '%')
                out1 = ''
                out2 = str(tmpread[0].split('\t')[4].strip()) + '\t' + str(chromocollect) + '\t' + \
                    str(actual_no_maps) + ' maps' + '\t' + ','.join(fullquerypercent) + '\t' + \
                    ','.join(evalue_total) + '\t' + ','.join(bitscore_ratio_total) + '\t' + ','.join(del_size) + ' ' + \
                    ','.join(del_range) + '\t' + ','.join(iso_ins_size) + ' ' + ','.join(iso_ins_range) + '\t' + \
                    ','.join(ins_size) + ' ' + ','.join(ins_range) + '\t' + ','.join(complex_ins_size) + '\t' + \
                    ','.join(tdup) + ' ' + ','.join(tdup_range) + '\t' + ','.join(inv) + ' ' + ','.join(inv_range) + '\t' + \
                    ','.join(intra_ins) + ' ' + ','.join(intra_ins_range) + '\t' + ','.join(inter_ins) + ' ' + \
                    ','.join(inter_ins_range) + '\t' + ','.join(intertx) + ' ' + ','.join(intertx_break) + '\t' + \
                    ','.join(strand_total) + '\t' + str(short) + '\t' + ','.join(complex_sv) + '\t' + \
                    ', '.join(total_sv_complex.split()) + '\t' + ','.join(sv_range) + '\t' + str(query) + '\t' + \
                    ','.join(piden_total) + '\t' + ','.join(mismatch_ratio_total) + '\t' + ','.join(gap_ratio_total)
            else:
                out1, out2 = '', ''
        else:
            out1, out2 = '', ''
    return out1, out2


# Out1 columns: 1=readname, 2=readlength, 3=total_chromosome(s), 4=No._of_maps, 5=query_signature, 6=evalue, 7=bitscore,
# 8=deletion, 9=Isolated_insertion, 10=Insertion, 11=complex_insertion, 12=tandemDupl, 13=inversion, 14=intra_insertion,
# 15=inter_insertion, 16=inter_tx, 17=strandness, 18=short_frag_perc, 19=satellite_perc, 20=complex_SV, 21=Total_SV_complex,
# 22=sv_range, 23=chromosome_order, 24=subject_string, 25=query_string, 26=%identity, 27=mismatch, 28=gap_ratio

# out1 = str(tmpread[0].split('\t')[4].strip()) + '\t' + str(readlength) + '\t' + str(chromocollect) + '\t' + \
#     str(actual_no_maps) + ' maps' + '\t' + ','.join(fullquerypercent) + '\t' + ','.join(evalue_total) + '\t' + \
#     ','.join(bitscore_ratio_total) + '\t' + ','.join(del_size) + ' ' + ','.join(del_range) + '\t' + \
#     ','.join(iso_ins_size) + ' ' + ','.join(iso_ins_range) + '\t' + ','.join(ins_size) + ' ' + \
#     ','.join(ins_range) + '\t' + ','.join(complex_ins_size) + '\t' + ','.join(tdup) + ' ' + \
#     ','.join(tdup_range) + '\t' + ','.join(inv) + ' ' + ','.join(inv_range) + '\t' + ','.join(intra_ins) + ' ' + \
#     ','.join(intra_ins_range) + '\t' + ','.join(inter_ins) + ' ' + ','.join(inter_ins_range) + '\t' + \
#     ','.join(intertx) + ' ' + ','.join(intertx_break) + '\t' + ','.join(strand_total) + '\t' + \
#     str(short) + '\t' + sat + '\t' + ','.join(complex_sv) + '\t' + ", ".join(total_sv_complex.split()) + '\t' + \
#     ','.join(sv_range) + '\t' + str(chromorder) + '\t' + ''.join(str(subject)) + '\t' + str(query) + '\t' + \
#     ','.join(piden_total) + '\t' + ','.join(mismatch_ratio_total) + '\t' + ','.join(gap_ratio_total)

# Out2: 0=readname, 1=total_chromosome(s), 2=No._of_maps, 3=query_signature, 4=evalue, 5=bitscore, 6=deletion,
# 7=Isolated_insertion, 8=Insertion, 9=complex_insertion, 10=tandemDupl, 11=inversion, 12=intra_insertion, 13=inter_insertion,
# 14=inter_tx, 15=strandness, 16=short_frag_perc, 17=complex_SV, 18=Total_SV_complex, 19=sv_range, 20=query_string, 21=%identity,
# 22=mismatch, 23=gap_ratio


# Check if single alignment covers (1-splitpct) of read length
def single_alignment_check(subdata, splitpct, minalign):
    align_percent = float(subdata[0].split('\t')[6]) / float(subdata[0].split('\t')[8])
    align_len = int(subdata[0].split('\t')[6])
    if align_percent >= (float(1)-splitpct) or align_len < minalign:
        return False
    else:
        return True


# Reduce mapping of multi-mapping reads
def reduce_multi_map(subdata):
    if len(subdata) > 4:
        sortdictlen = {}
        for line in subdata:
            sortdictlen[line] = int(line.split('\t')[6])
        # Sorting according to alignment length
        subdata2 = [key for (key, value) in sorted(sortdictlen.items(), key=lambda x: x[1], reverse=True)]
        subdata3 = subdata2[0:4]  # Selecting top 4 longest alignments
        # Sort back according to query start
        sortdict = {}
        for e in subdata3:
            sortdict[int(e.split('\t')[5])] = e
        subdata4 = [value for (key, value) in sorted(sortdict.items())]
        return subdata4
    else:
        return subdata


# For future implementation to limit multi-mapping
def discard_multi_map(actual_no_maps, multi_limit):
    # Throw away multi-mapper reads
    if int(actual_no_maps) >= multi_limit:
        return False
    else:
        return True

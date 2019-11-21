"""
Functions for SV breakend clustering.

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


from collections import OrderedDict
from pybedtools import BedTool


def sv_cluster(subdata, parse, buf, maxovl):
    rangedict, svsizedictsort, infodict, bed_str1, bed_str2 = rangecollect(parse, buf)
    bed1 = BedTool(bed_str1, from_string=True)
    bed1 = bed1.sort()
    bed2 = BedTool(bed_str2, from_string=True)
    bed2 = bed2.sort()
    insect = bed1.intersect(bed2, wa=True, wb=True)
    intersect = []
    for line in insect:
        line = list(line)
        if line[3] != line[8]:
            intersect.append('\t'.join(line))
    connectdict, classdict = intersection(intersect, bed1)
    tier1, tier2, tier3 = classranker(classdict, svsizedictsort)
    readdicta, readdictb = connector(connectdict, tier1, tier2, tier3, classdict)
    readdictc = classupdate(readdictb)
    newdict = bestread(readdicta, classdict, readdictc)
    bed_str3 = normalbed(subdata)
    bed_str4 = svbed(newdict, infodict)
    bed3 = BedTool(bed_str3, from_string=True)
    bed3 = bed3.sort()
    bed4 = BedTool(bed_str4, from_string=True)
    bed4 = bed4.sort()
    insect2 = bed3.intersect(bed4, wa=True, wb=True)
    intersect2 = []
    for line in insect2:
        line = list(line)
        if line[3] != line[7]:
            intersect2.append('\t'.join(line))
    newdictnouniq = remove_uniq(newdict)  # Remove unique identifier from read name
    svnormalcov = intersection2(intersect2, newdict, newdictnouniq, infodict, classdict)
    output = arrange(newdict, infodict, svnormalcov, maxovl)
    return output


# Function to convert non-digits to zero
def alpha(n):
    if str(n).isdigit():
        return int(n)
    else:
        return int(0)


# Function to convert . to zero
def dotter(dot):
    if dot == '.':
        return float(0)
    else:
        return float(dot)


# Function to set minimum value to 0
def mz(value):
    if value < 0:
        return 0
    else:
        return value


# Function to collect information of each breakpoint entry (e.g. chromosome, range, SV size, SV type)
def rangecollect(x, buf):
    namerepeat = OrderedDict()
    rangedict = OrderedDict()
    svsizedict = OrderedDict()
    infodict = OrderedDict()
    bed1 = []
    bed2 = []
    for i in x:
        try:
            infodict[i.split('\t')[8]].append('\t'.join(i.split('\t')[0:10]))
        except KeyError:
            infodict[i.split('\t')[8]] = []
            infodict[i.split('\t')[8]].append('\t'.join(i.split('\t')[0:10]))
        try:
            namerepeat[i.split('\t')[8]] == 1
        except KeyError:
            namerepeat[i.split('\t')[8]] = 1
            rnameidx = i.split('\t')[8]
            svtype = str(svtypecorrector(i.split('\t')[3].split(' ')[0]))
            if len(i.split('\t')[6].split('~')) == 2:  # Not Inter translocation
                if len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 1:
                    chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                    le = int(i.split('\t')[6].split('~')[1].split(':')[1])
                    svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[-1])
                    bed1.append(chm1 + '\t' + str(le) + '\t' + str(le+1) + '\t' + rnameidx + '-l' + '\t' + svtype)
                    bed2.append(chm1 + '\t' + str(mz(le-buf)) + '\t' + str(le+buf) + '\t' + rnameidx + '-l' + '\t' + svtype)
                elif len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 2:
                    chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                    le = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
                    r = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
                    svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[-1])
                    bed1.append(chm1 + '\t' + str(le) + '\t' + str(le+1) + '\t' + rnameidx + '-l' + '\t' + svtype)
                    bed1.append(chm1 + '\t' + str(r) + '\t' + str(r+1) + '\t' + rnameidx + '-r' + '\t' + svtype)
                    bed2.append(chm1 + '\t' + str(mz(le-buf)) + '\t' + str(le+buf) + '\t' + rnameidx + '-l' + '\t' + svtype)
                    bed2.append(chm1 + '\t' + str(mz(r-buf)) + '\t' + str(r+buf) + '\t' + rnameidx + '-r' + '\t' + svtype)
            elif len(i.split('\t')[6].split('~')) == 3:  # Inter translocation
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                chm2 = i.split('\t')[6].split('~')[2].split(':')[0]
                le = int(i.split('\t')[6].split('~')[1].split(':')[1])
                r = int(i.split('\t')[6].split('~')[2].split(':')[1])
                svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[-1])
                bed1.append(chm1 + '\t' + str(le) + '\t' + str(le+1) + '\t' + rnameidx + '-l' + '\t' + svtype)
                bed1.append(chm2 + '\t' + str(r) + '\t' + str(r+1) + '\t' + rnameidx + '-r' + '\t' + svtype)
                bed2.append(chm1 + '\t' + str(mz(le-buf)) + '\t' + str(le+buf) + '\t' + rnameidx + '-l' + '\t' + svtype)
                bed2.append(chm2 + '\t' + str(mz(r-buf)) + '\t' + str(r+buf) + '\t' + rnameidx + '-r' + '\t' + svtype)
    svsizedictsort = [key for (key, value) in sorted(svsizedict.items(), key=lambda y: y[1], reverse=True)]
    bed_str1 = '\n'.join(bed1)
    bed_str2 = '\n'.join(bed2)
    return rangedict, svsizedictsort, infodict, bed_str1, bed_str2


# Function to standardize names of SV types
def svtypecorrector(svtype):
    if svtype == 'S-Nov_Ins_bp':
        return 'bp_Nov_Ins'
    elif svtype == 'E-Nov_Ins_bp':
        return 'bp_Nov_Ins'
    elif svtype == 'Inter-Ins(1)':
        return 'Inter-Ins'
    elif svtype == 'Inter-Ins(2)':
        return 'Inter-Ins'
    elif svtype == 'InterTx':
        return 'Inter'
    elif svtype == 'Inv(1)':
        return 'Inv2'
    elif svtype == 'Inv(2)':
        return 'Inv2'
    elif svtype == 'Intra-Ins(1)':
        return 'Intra-Ins2'
    elif svtype == 'Intra-Ins(2)':
        return 'Intra-Ins2'
    else:
        return svtype


def svgroup(svclass):
    if svclass == 'bp_Nov_Ins':
        return ['bp_Nov_Ins', 'Nov_Ins', 'Inter-Ins', 'Inter', 'Inv2', 'Inv', 'Intra-Ins2', 'Intra-Ins', 'Del', 'TDupl']
    elif svclass == 'Nov_Ins':
        return ['bp_Nov_Ins', 'Nov_Ins']
    elif svclass == 'Inter-Ins':
        return ['bp_Nov_Ins', 'Inter-Ins', 'Inter']
    elif svclass == 'Inter':
        return ['bp_Nov_Ins', 'Inter-Ins', 'Inter']
    elif svclass == 'Inv2':
        return ['bp_Nov_Ins', 'Inv2', 'Inv']
    elif svclass == 'Inv':
        return ['bp_Nov_Ins', 'Inv2', 'Inv']
    elif svclass == 'Intra-Ins2':
        return ['bp_Nov_Ins', 'Intra-Ins2', 'Intra-Ins']
    elif svclass == 'Intra-Ins':
        return ['bp_Nov_Ins', 'Intra-Ins2', 'Intra-Ins']
    elif svclass == 'Del':
        return ['bp_Nov_Ins', 'Del']
    elif svclass == 'TDupl':
        return ['bp_Nov_Ins', 'TDupl']


# Function to connect overlaping overlaping breakpoint entries
def intersection(intersect, bed1):
    connectdict = OrderedDict()
    classdict = OrderedDict()
    for i in bed1:
        i = list(i)
        connectdict[i[3]] = []
        classdict[i[3][0:-2]] = i[4]
    for i in intersect:
        if i.split('\t')[4] in svgroup(i.split('\t')[9]):
            connectdict[i.split('\t')[8]].append(i.split('\t')[3][0:-2])
    return connectdict, classdict


# Remove unique identifier from read name
def remove_uniq(x):
    newdictnouniq = OrderedDict()
    for i in x:
        newdictnouniq[i] = []
        for j in x[i]:
            newdictnouniq[i].append(j[:-8])
    return newdictnouniq


def intersection2(intersect2, x, newdictnouniq, infodict, classdict):
    svnormalconnect = OrderedDict()
    svnormalcov = OrderedDict()
    for key in x:
        svnormalconnect[key] = []
    for i in intersect2:
        if i.split('\t')[3] not in newdictnouniq[i.split('\t')[7]]:
            if classdict[i.split('\t')[7]] == 'TDupl':  # If Dup, make sure normal read covers both breakpoints of Dup
                # for n in infodict[i.split('\t')[7]]:  # Redundant as each novel adjacency on read can have max 2 breakpoints
                if int(i.split('\t')[1]) <= int(infodict[i.split('\t')[7]][0].split('\t')[5]) <= int(i.split('\t')[2]) and \
                        int(i.split('\t')[1]) <= int(infodict[i.split('\t')[7]][1].split('\t')[5]) <= int(i.split('\t')[2]):
                    svnormalconnect[i.split('\t')[7]].append(i.split('\t')[3])
            else:
                svnormalconnect[i.split('\t')[7]].append(i.split('\t')[3])
    for key in svnormalconnect:
        svnormalcov[key] = float(len(set(svnormalconnect[key])))  # set() removes duplicates
    return svnormalcov


# Ordering classes
def classranker(x, svsizedictsort):
    tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
    tier2sv = ['Inv', 'Del', 'TDupl', 'Nov_Ins', 'Intra-Ins', 'Inter']
    tier1sv = ['bp_Nov_Ins']
    tier1 = []
    tier2 = []
    tier3 = []
    for keys in svsizedictsort:
        if x[keys] in tier1sv:
            tier1.append(keys)
        elif x[keys] in tier2sv:
            tier2.append(keys)
        elif x[keys] in tier3sv:
            tier3.append(keys)
        else:
            raise Exception("Error: SV class error")
    return tier1, tier2, tier3


def connector(x, tier1, tier2, tier3, classdict):
    namerepeat = OrderedDict()
    readdicta = OrderedDict()
    readdictb = OrderedDict()
    alltier = tier3 + tier2 + tier1
    for keys in alltier:
        if keys not in namerepeat:
            namerepeat[keys] = 1
            readdicta[keys] = []
            readdictb[keys] = []
            readdictb[keys].append(classdict[keys])
            for i in x[keys + '-l']:
                if i not in namerepeat:
                    if keys + '-r' in x:
                        if i in x[keys + '-r']:
                            namerepeat[i] = 1
                            readdicta[keys].append(i)
                            readdictb[keys].append(classdict[i])
                        else:
                            if i in tier1:  # Only for single-bp-end reads
                                namerepeat[i] = 1
                                readdicta[keys].append(i)
                                readdictb[keys].append(classdict[i])
                    else:
                        namerepeat[i] = 1
                        readdicta[keys].append(i)
                        readdictb[keys].append(classdict[i])
            if keys + '-r' in x:
                for i in x[keys + '-r']:
                    if i not in namerepeat:
                        if i in tier1:  # Only for single-bp-end reads
                            namerepeat[i] = 1
                            readdicta[keys].append(i)
                            readdictb[keys].append(classdict[i])
            for i in readdicta[keys]:
                for k in x[i + '-l']:
                    if k not in namerepeat:
                        if i + '-r' in x:
                            if k in x[i + '-r']:
                                namerepeat[k] = 1
                                readdicta[keys].append(k)
                                readdictb[keys].append(classdict[k])
                        else:
                            pass  # Ignoring fellows from single-bp-end reads
    return readdicta, readdictb


# Function to rank and filter best sv type
def classupdate(readdictb):
    readdictc = OrderedDict()
    d = OrderedDict()
    tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
    tier2sv = ['Inv', 'Del', 'TDupl', 'Nov_Ins', 'Intra-Ins', 'Inter']
    tier1sv = ['bp_Nov_Ins']
    for keys in readdictb:
        for i in readdictb[keys]:
            if i not in d:
                d[i] = 1
            else:
                d[i] += 1
        if any(y in tier3sv for y in d):
            for s in tier2sv:
                try:
                    del d[s]
                except KeyError:
                    pass
            try:
                del d["bp_Nov_Ins"]
            except KeyError:
                pass
        elif any(y in tier2sv for y in d):
            try:
                del d["bp_Nov_Ins"]
            except KeyError:
                pass
        elif any(y in tier1sv for y in d):
            pass
        else:
            raise Exception("Error: Unknown SV class")
        mainsv = [key for (key, value) in sorted(d.items(), key=lambda x:x[1], reverse=True)][0]
        readdictc[keys] = mainsv
        d = OrderedDict()
    return readdictc


# Function to find the lead breakpoint entry based on main sv type first come across
def bestread(x, classdict, readdictc):
    newdict = OrderedDict()
    for keys in x:
        main = ''
        mainsv = readdictc[keys]
        if classdict[keys] == mainsv:
            main = keys
            newdict[main] = []
            for k in x[main]:
                newdict[main].append(k)
        else:
            for k in x[keys]:
                if classdict[k] == mainsv:
                    main = k
                    break
            newdict[main] = []
            newdict[main].append(keys)
            for r in x[keys]:
                if r != main:
                    newdict[main].append(r)
    return newdict


# Function to generate bed file from hsblast subdata
def normalbed(subdata):
    totalbed = []
    len_buf = 400  # Arbituary length deduction buffer for normal reads, increases stringency
    for line in subdata:
        coord1 = int(line.split('\t')[1]) + len_buf
        coord2 = int(line.split('\t')[1]) + int(line.split('\t')[2]) - len_buf
        if coord2 - coord1 > 0:
            totalbed.append(line.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + line.split('\t')[4])
        else:
            pass
    bed_str3 = '\n'.join(totalbed)
    return bed_str3


def svbed(x, infodict):
    totalsvbed = []
    for key in x:
        for n in infodict[key]:
            totalsvbed.append('\t'.join(n.split('\t')[4:6]) + '\t' + str(int(n.split('\t')[5]) + 1) + '\t' + key)
    bed_str4 = '\n'.join(totalsvbed)
    return bed_str4


# Function to count total number of unique reads in an overlaping breakpoint entry
def countcov(key, x):
    readdict = OrderedDict()
    readdict[key[:-6]] = 1
    for i in x[key]:
        readdict[i[:-6]] = 1
    return len(readdict)


# Function to parse lead overlaping breakpoint entry and long read coverage into output
def arrange(x, infodict, svnormalcov, maxovl):
    output = []
    for key in x:
        lcov = countcov(key, x)
        if lcov == 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t.\t' + str(svnormalcov[key]))
        elif 1 < lcov < maxovl:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t' + ','.join(x[key]) + '\t' + str(svnormalcov[key]))
    return output


# Function to parse lead overlaping breakpoint entry, long read coverage and average DNN score into output (Not in use)
def arrangednn(x, infodict, svnormalcov, dnndict):
    output = []
    for key in x:
        dnn = avgdnn(key, x[key], dnndict)
        lcov = countcov(key, x)
        if lcov == 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t.\t' + str(svnormalcov[key]) + '\t' + str(dnn) + '\n')
        elif lcov > 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t' + ','.join(x[key]) + '\t' + str(svnormalcov[key]) + '\t' + str(dnn) +
                              '\n')
    return output


# Function to calculate the average short read coverage of an overlaping breakpoint entry (Not in use)
def avgdnn(key, value, dnndict):
    dnn = 0
    num = len(value) + 1
    dnn = dnn + dnndict[key]
    for k in value:
        dnn = dnn + dnndict[k]
    avg = float(dnn)/num
    return avg


# Note: Some bps seen in parse_file would be missing
# in overlap_file because they are overlapped by another
# bp within the same read

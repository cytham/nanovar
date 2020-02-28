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


# SV clustering
def sv_cluster(subdata, parse, buf, maxovl, mincov, ins_switch):
    readteam, infodict, classdict, mainclass, svsizedict = rangecollect(parse, buf, ins_switch)
    bed_str3 = normalbed(subdata)
    bed_str4 = svbed(readteam, mainclass)
    bed3 = BedTool(bed_str3, from_string=True)
    bed3 = bed3.sort()
    bed4 = BedTool(bed_str4, from_string=True)
    bed4 = bed4.sort()
    insect2 = bed3.intersect(bed4, wa=True, wb=True)
    intersect2 = []
    for line in insect2:
        line = list(line)
        intersect2.append('\t'.join(line))
    newdictnouniq = remove_uniq(readteam)  # Remove unique identifier from read name
    svnormalcov = intersection2(intersect2, readteam, newdictnouniq, mainclass)
    output = arrange(svnormalcov, readteam, maxovl, mincov, infodict, mainclass, svsizedict)
    return output


# Function to convert non-digits to zero
def alpha(n):
    if str(n).isdigit():
        return int(n)
    else:
        return int(0)


# Function to collect information of each breakpoint entry and cluster breakpoints
def rangecollect(x, buf, ins_switch):
    classdict = OrderedDict()
    svsizedict = OrderedDict()
    infodict = OrderedDict()
    classclust = OrderedDict()
    for i in x:
        rnameidx = i.split('\t')[8]
        infodict[rnameidx] = '\t'.join(i.split('\t')[0:10])
        svtype = str(svtypecorrector(i.split('\t')[3].split(' ')[0]))
        classdict[rnameidx] = svtype
        svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[1].split('~')[0])
        # if svtype != 'bp_Nov_Ins':
        if len(i.split('\t')[6].split('~')) == 2:  # Not Inter translocation
            if len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 1:
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                le = int(i.split('\t')[6].split('~')[1].split(':')[1])
                if chm1 not in classclust:
                    classclust[chm1] = OrderedDict()
                classclust[chm1][rnameidx + '-l'] = le
            elif len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 2:
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                le = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
                r = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
                if chm1 not in classclust:
                    classclust[chm1] = OrderedDict()
                classclust[chm1][rnameidx + '-l'] = le
                classclust[chm1][rnameidx + '-r'] = r
        elif len(i.split('\t')[6].split('~')) == 3:  # Inter translocation
            chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
            chm2 = i.split('\t')[6].split('~')[2].split(':')[0]
            le = int(i.split('\t')[6].split('~')[1].split(':')[1])
            r = int(i.split('\t')[6].split('~')[2].split(':')[1])
            if chm1 not in classclust:
                classclust[chm1] = OrderedDict()
            classclust[chm1][rnameidx + '-l'] = le
            if chm2 not in classclust:
                classclust[chm2] = OrderedDict()
            classclust[chm2][rnameidx + '-r'] = r
    chrnamelist = OrderedDict()
    chrcoorlist = OrderedDict()
    for chrm in classclust:
        classclust[chrm] = OrderedDict(sorted(classclust[chrm].items(), key=lambda y: y[1]))
        chrnamelist[chrm] = []
        chrcoorlist[chrm] = []
        for key in classclust[chrm]:
            chrnamelist[chrm].append(key[:-2])
            chrcoorlist[chrm].append(classclust[chrm][key])
    readteam, mainclass = cluster(chrnamelist, chrcoorlist, buf, svsizedict, classdict, ins_switch)
    return readteam, infodict, classdict, mainclass, svsizedict


# Cluster coordinates
def cluster(chrnamelist, chrcoorlist, buf, svsizedict, classdict, ins_switch):
    cluster2read = OrderedDict()
    read2cluster = OrderedDict()
    for chrm in chrnamelist:
        n = len(chrnamelist[chrm])
        last = None
        group = []
        groupnames = []
        for i in range(n):
            if last is None or chrcoorlist[chrm][i] - last <= buf:
                group.append(chrcoorlist[chrm][i])
                groupnames.append(chrnamelist[chrm][i])
            else:
                avgcoord = int(round(sum(group) / len(group), 0))
                cluster2read[chrm + ':' + str(avgcoord)] = sorted(set(groupnames))
                for read in sorted(set(groupnames)):
                    if read not in read2cluster:
                        read2cluster[read] = []
                    read2cluster[read].append(chrm + ':' + str(avgcoord))
                group = [chrcoorlist[chrm][i]]
                groupnames = [chrnamelist[chrm][i]]
            last = chrcoorlist[chrm][i]
        if group:
            avgcoord = int(round(sum(group) / len(group), 0))
            cluster2read[chrm + ':' + str(avgcoord)] = sorted(set(groupnames))
            for read in sorted(set(groupnames)):
                if read not in read2cluster:
                    read2cluster[read] = []
                read2cluster[read].append(chrm + ':' + str(avgcoord))
    clusterlink = OrderedDict()
    coordlist2 = []
    for read in read2cluster:
        if read2cluster[read][0] not in clusterlink:  # priming
            clusterlink[read2cluster[read][0]] = []
        if len(read2cluster[read]) == 1:
            pass
        elif len(read2cluster[read]) == 2:
            clusterlink[read2cluster[read][0]].append(read2cluster[read][1])
            coordlist2.append(read2cluster[read][1])
        else:
            raise Exception("Read %s has none or more than 2 clusters" % read)
    for clusters in coordlist2:  # for clusters that have been right paired
        try:
            if not clusterlink[clusters]:  # if that cluster has no right paired/is single cluster
                del clusterlink[clusters]  # delete the cluster
        except KeyError:
            pass
    for clusters in clusterlink:
        clusterlink[clusters] = sorted(set(clusterlink[clusters]))
    readteam = OrderedDict()
    mainclass = OrderedDict()
    pastbps = {}
    for coord in cluster2read:  # For each clusters
        if coord in clusterlink:  # If cluster is either unpaired or paired (left)
            if not clusterlink[coord]:  # If cluster is unpaired, svclass concordance matters here
                leaders, teams, svclasses = svclass_sep(cluster2read[coord], classdict, svsizedict)
                # lead, main = leadread(cluster2read[coord], svsizedict, classdict, ins_switch)
                no = len(leaders)
                for i in range(no):
                    if svclasses[i] in ['Nov_Ins', 'bp_Nov_Ins']:
                        readteam[':'.join([coord.split(':')[0], str(int(coord.split(':')[1])+i)])] = [leaders[i]]
                        readteam[':'.join([coord.split(':')[0], str(int(coord.split(':')[1])+i)])].extend(sorted(set(
                            teams[i]).difference([leaders[i]])))
                        mainclass[':'.join([coord.split(':')[0], str(int(coord.split(':')[1])+i)])] = svclasses[i]
                    else:
                        readteam[':'.join([coord.split(':')[0], str(int(coord.split(':')[1])+i)]) + '-' + ':'.join([
                            coord.split(':')[0], str(int(coord.split(':')[1])+1+i)])] = [leaders[i]]
                        readteam[':'.join([coord.split(':')[0], str(int(coord.split(':')[1])+i)]) + '-' + ':'.join([
                            coord.split(':')[0], str(int(coord.split(':')[1])+1+i)])].extend(
                            sorted(set(teams[i]).difference([leaders[i]])))
                        mainclass[':'.join([coord.split(':')[0], str(int(coord.split(':')[1])+i)]) + '-' + ':'.join([
                            coord.split(':')[0], str(int(coord.split(':')[1])+1+i)])] = svclasses[i]
            else:  # If cluster is paired (left)
                for coord2 in clusterlink[coord]:
                    readlist = []
                    readlist.extend(sorted(set(cluster2read[coord]).intersection(cluster2read[coord2])))
                    for read in sorted(set(cluster2read[coord]).difference(cluster2read[coord2])):
                        if classdict[read] == 'bp_Nov_Ins':
                            if read not in pastbps:
                                readlist.append(read)
                                pastbps[read] = 1
                    lead, main = leadread(sorted(set(readlist)), svsizedict, classdict, ins_switch)
                    readteam[coord + '-' + coord2] = [lead]
                    readteam[coord + '-' + coord2].extend(sorted(set(readlist).difference([lead])))
                    mainclass[coord + '-' + coord2] = main
    return readteam, mainclass


def leadread(reads, svsizedict, classdict, ins_switch):
    if ins_switch:
        mainsvclass = mainclasssv_ins(reads, classdict)
    else:
        mainsvclass = mainclasssv(reads, classdict)
    sizedict = OrderedDict()
    leader = ''
    for read in reads:
        sizedict[read] = svsizedict[read]
    sizedictsort = [key for (key, value) in sorted(sizedict.items(), key=lambda y: y[1], reverse=True)]
    for read in sizedictsort:
        if classdict[read] == mainsvclass:
            leader = read
            break
    if leader == '':
        raise Exception("Error: Main SV class not found")
    return leader, mainsvclass


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


# Remove unique identifier from read name
def remove_uniq(readteam):
    newdictnouniq = OrderedDict()
    for clusters in readteam:
        newdictnouniq[clusters] = []
        for read in readteam[clusters]:
            newdictnouniq[clusters].append(read[:-6])
    return newdictnouniq


# Function to rank and filter best sv type
def mainclasssv(reads, classdict):
    tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
    tier2sv = ['Inv', 'Del', 'TDupl', 'Nov_Ins', 'Intra-Ins', 'Inter']
    tier1sv = ['bp_Nov_Ins']
    tmpdict = OrderedDict()
    for read in reads:
        svclass = classdict[read]
        if svclass not in tmpdict:
            tmpdict[svclass] = 1
        else:
            tmpdict[svclass] += 1
    if any(y in tier3sv for y in tmpdict):
        for s in tier2sv:
            try:
                del tmpdict[s]
            except KeyError:
                pass
        try:
            del tmpdict["bp_Nov_Ins"]
        except KeyError:
            pass
    elif any(y in tier2sv for y in tmpdict):
        try:
            del tmpdict["bp_Nov_Ins"]
        except KeyError:
            pass
    elif any(y in tier1sv for y in tmpdict):
        pass
    else:
        raise Exception("Error: Unknown SV class")
    mainsvclass = [key for (key, value) in sorted(tmpdict.items(), key=lambda x: x[1], reverse=True)][0]
    return mainsvclass


# Function to rank and filter SV type after INS re-analysis
def mainclasssv_ins(reads, classdict):
    tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2', 'TDupl']
    tier2sv = ['Inv', 'Del', 'Nov_Ins', 'Intra-Ins', 'Inter']
    tier1sv = ['bp_Nov_Ins']
    tmpdict = OrderedDict()
    for read in reads:
        svclass = classdict[read]
        if svclass not in tmpdict:
            tmpdict[svclass] = 1
        else:
            tmpdict[svclass] += 1
    if any(y in tier3sv for y in tmpdict):
        for s in tier2sv:
            try:
                del tmpdict[s]
            except KeyError:
                pass
        try:
            del tmpdict["bp_Nov_Ins"]
        except KeyError:
            pass
    elif any(y in tier2sv for y in tmpdict):
        try:
            del tmpdict["bp_Nov_Ins"]
        except KeyError:
            pass
    elif any(y in tier1sv for y in tmpdict):
        pass
    else:
        raise Exception("Error: Unknown SV class")
    mainsvclass = [key for (key, value) in sorted(tmpdict.items(), key=lambda x: x[1], reverse=True)][0]
    return mainsvclass


# Segregate an unpaired cluster by different sv classes
def svclass_sep(reads, classdict, svsizedict):
    svreadteam = OrderedDict()
    for read in reads:
        svclass = classdict[read]
        if svclass == 'bp_Nov_Ins':  # Bundle bp_Nov_Ins and Nov_Ins together
            svclass = 'Nov_Ins'
        try:
            svreadteam[svclass].append(read)
        except KeyError:
            svreadteam[svclass] = [read]
    leaders = []
    teams = []
    svclasses = []
    for svtype in svreadteam:
        sizedict = OrderedDict()
        leader = ''
        for read in svreadteam[svtype]:
            sizedict[read] = svsizedict[read]
        sizedictsort = [key for (key, value) in sorted(sizedict.items(), key=lambda y: y[1], reverse=True)]
        for read in sizedictsort:
            if classdict[read] == svtype:
                leader = read
                break
        if leader == '':
            for read in sizedictsort:
                if classdict[read] == 'bp_Nov_Ins':
                    leader = read
                    break
            if leader == '':
                raise Exception("Error: Main SV seg class not found")
        leaders.append(leader)
        teams.append(svreadteam[svtype])
        svclasses.append(classdict[leader])
    return leaders, teams, svclasses


# Convert single element list to .
def dotter(readlist):
    if len(readlist) == 1:
        return '.'
    else:
        return readlist[1:]


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


# Function to generate bed file for SV breakpoints
def svbed(readteam, mainclass):
    totalsvbed = []
    for clusters in readteam:
        svtype = mainclass[clusters]
        clusterlist = clusters.split('-')
        for n in clusterlist:
            totalsvbed.append(n.split(':')[0] + '\t' + n.split(':')[1] + '\t' + str(int(n.split(':')[1]) + 1) + '\t' +
                              clusters + '\t' + svtype)
    bed_str4 = '\n'.join(totalsvbed)
    return bed_str4


# Parse BED intersection output
def intersection2(intersect2, readteam, newdictnouniq, mainclass):
    svnormalconnect = OrderedDict()
    svnormalcov = OrderedDict()
    for clusters in readteam:
        svtype = mainclass[clusters]
        if svtype in ['Intra-Ins2', 'Intra-Ins', 'Inter', 'Inter-Ins']:
            svnormalconnect[clusters] = [[], []]
        else:
            svnormalconnect[clusters] = []
    for line in intersect2:
        clusters = line.split('\t')[7]
        svtype = line.split('\t')[8]
        # hybrid = clusters + '~' + svtype
        read = line.split('\t')[3]
        coord1 = int(line.split('\t')[1])
        coord2 = int(line.split('\t')[2])
        if read not in newdictnouniq[clusters]:
            if svtype == 'TDupl':  # If Dup, make sure normal read covers both breakpoints of Dup
                if coord1 <= int(clusters.split('-')[0].split(':')[1]) <= coord2 and \
                        coord1 <= int(clusters.split('-')[1].split(':')[1]) <= coord2:
                    svnormalconnect[clusters].append(read)
            elif svtype in ['Intra-Ins2', 'Intra-Ins', 'Inter', 'Inter-Ins']:
                if line.split('\t')[5] == clusters.split('-')[0].split(':')[1]:  # Figuring which breakend (left)
                    svnormalconnect[clusters][0].append(read)
                elif line.split('\t')[5] == clusters.split('-')[1].split(':')[1]:  # Figuring which breakend (right)
                    svnormalconnect[clusters][1].append(read)
                else:
                    raise Exception('Error: BND coordinate not in records')
            else:
                svnormalconnect[clusters].append(read)
    for clusters in svnormalconnect:
        svtype = mainclass[clusters]
        if svtype in ['Del', 'Inv']:
            svnormalcov[clusters] = int(round(len(svnormalconnect[clusters])/2, 0))
        elif svtype in ['Intra-Ins2', 'Intra-Ins', 'Inter', 'Inter-Ins']:
            svnormalcov[clusters] = min(len(set(svnormalconnect[clusters][0])), len(set(svnormalconnect[clusters][1])))
        else:
            svnormalcov[clusters] = int(len(set(svnormalconnect[clusters])))  # set() removes duplicates
    return svnormalcov


# Function to count total number of unique reads in an overlaping breakpoint entry
def countcov(readlist):
    newlist = []
    for read in readlist:
        newlist.append(read[:-6])
    return len(set(newlist))


# Parse info into output
def arrange(svnormalcov, readteam, maxovl, mincov, infodict, mainclass, svsizedict):
    output = []
    n = 1
    for clusters in readteam:
        svtype = mainclass[clusters]
        lcov = countcov(readteam[clusters])
        bestread = readteam[clusters][0]
        if mincov <= lcov <= maxovl and svnormalcov[clusters] <= maxovl:
            if svtype in ['Inter', 'Inter-Ins']:
                output.append(
                    '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                    infodict[bestread].split('\t')[6].split('~')[0] +
                    '~' + clusters.split('-')[0] + '~' + clusters.split('-')[1] + '\t' +
                    '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                    ','.join(dotter(readteam[clusters])) + '\t' + str(svnormalcov[clusters])
                )
            elif svtype in ['Nov_Ins', 'Del'] and svsizedict[bestread] <= 200 and lcov < max(mincov, 2):
                continue
            elif svtype == 'bp_Nov_Ins' and lcov < max(mincov, 2):
                continue
            else:
                if len(clusters.split('-')) == 1:
                    output.append(
                        '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                        infodict[bestread].split('\t')[6].split('~')[0] + '~' + clusters + '-' +
                        str(int(clusters.split(':')[1]) + 1) + '\t' +
                        '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                        ','.join(dotter(readteam[clusters])) + '\t' + str(svnormalcov[clusters])
                    )
                elif len(clusters.split('-')) == 2:
                    output.append(
                        '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                        infodict[bestread].split('\t')[6].split('~')[0] + '~' +
                        clusters.split('-')[0] + '-' + clusters.split('-')[1].split(':')[1] + '\t' +
                        '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                        ','.join(dotter(readteam[clusters])) + '\t' + str(svnormalcov[clusters])
                    )
                else:
                    raise Exception('Error: Cluster name error')
            n += 1
    return output


# Note: Some bps seen in parse_file would be missing
# in cluster_file because they are overlapped by another
# bp within the same read

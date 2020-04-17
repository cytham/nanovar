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
from collections import OrderedDict, defaultdict
from pybedtools import BedTool


# SV clustering main function
def sv_cluster(
        subdata,
        parse,
        buf,
        maxovl,
        mincov,
        contigs,
        ins_switch
):
    """Returns a list of clustered SV breakends in strings

Keyword arguments:
- subdata: BED information for all reads in input
- parse: parse data from breakpoint_parser
- buf: coordinate clustering buffer in base pairs
- maxovl: maximum number of read overlap allowed for breakend-supporting or opposing- reads
- mincov: minimum number of breakend-supporting reads required for a breakend
- contigs: list of contig names from input reference genome
- ins_switch: boolean to switch between pre and post blast re-evaluation

    """
    # Parse SV locations and cluster them
    readteam, infodict, classdict, mainclass, svsizedict = rangecollect(parse,
                                                                        buf,
                                                                        contigs,
                                                                        ins_switch,
                                                                        mincov,
                                                                        )

    # Convert subdata to BED format
    bed3 = normalbed(subdata)
    bed3 = BedTool('\n'.join(bed3), from_string=True)
    bed3 = bed3.sort()

    # Convert SV breakends to BED format
    bed4 = svbed2(readteam, mainclass)
    bed4 = BedTool('\n'.join(bed4), from_string=True)
    bed4 = bed4.sort()

    # Intersect all read regions to SV breakends
    intersect = bed3.intersect(bed4, wa=True, wb=True)

    # Create a duplicated readteam dict where reads have their unique identifier removed
    newdictnouniq = remove_uniq(readteam)

    # Parse intersection BED to discover breakend-opposing reads for each breakend
    svnormalcov = intersection3(intersect, readteam, newdictnouniq, mainclass)

    # Output all information into a list of strings
    output = arrange(svnormalcov, readteam, maxovl, mincov, infodict, mainclass, svsizedict)
    return output


# Function to collect information of each breakpoint entry and cluster breakpoints
def rangecollect(x, buf, contigs, ins_switch, mincov):

    classdict, svsizedict, infodict = OrderedDict(), OrderedDict(), OrderedDict()
    leftclust, rightclust = OrderedDict(), OrderedDict()

    for c in contigs:  # Might replace with defaultdict in the future
        leftclust[c] = OrderedDict()
        rightclust[c] = OrderedDict()

    for i in x:
        rnameidx = i.split('\t')[8]  # Read unique id
        infodict[rnameidx] = '\t'.join(i.split('\t')[0:10])  # Read info storage from total_out file
        svtype = str(svtypecorrector(i.split('\t')[3].split(' ')[0]))  # Get sv type and slightly modify it
        classdict[rnameidx] = svtype  # Assign sv type to classdict
        svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[1].split('~')[0])  # Assign sv size to svsizedict
        if len(i.split('\t')[6].split('~')) == 2:  # Not Inter translocation
            if svtype in ['Nov_Ins', 'bp_Nov_Ins']:  # If single breakend
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]  # Get contig name
                le = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])  # Take first coord of Nov_Ins or bp_Nov_Ins
                col4 = i.split('\t')[3]  # Get SV type label
                # Assign coord to left or right group depending on strand
                if col4.split(' ')[0][0] == 'S':  # Head novel sequence
                    if col4.split(' ')[1].split('~')[1] == '+':
                        rightclust[chm1][rnameidx + '-l'] = le
                    elif col4.split(' ')[1].split('~')[1] == '-':
                        leftclust[chm1][rnameidx + '-l'] = le
                elif col4.split(' ')[0][0] == 'E':  # Tail novel sequence
                    if col4.split(' ')[1].split('~')[1] == '+':
                        leftclust[chm1][rnameidx + '-l'] = le
                    elif col4.split(' ')[1].split('~')[1] == '-':
                        rightclust[chm1][rnameidx + '-l'] = le
                elif col4.split(' ')[0][0] == 'N':  # Nov insertion with both breakends of insertion event
                    # Assign single coord to both left and right groups
                    leftclust[chm1][rnameidx + '-l'] = le
                    rightclust[chm1][rnameidx + '-r'] = le
                else:
                    raise Exception('Error, single SV breakend has unknown naming')
            elif len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 2:  # If double breakends
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]  # Get contig name
                le = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])  # Get default left coord
                r = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[1])  # Get default right coord
                if svtype == 'TDupl':  # If tandem dup, reverse left and right groups
                    leftclust[chm1][rnameidx + '-l'] = r
                    rightclust[chm1][rnameidx + '-r'] = le
                else:
                    leftclust[chm1][rnameidx + '-l'] = le
                    rightclust[chm1][rnameidx + '-r'] = r
        elif len(i.split('\t')[6].split('~')) == 3:  # Inter translocation
            chm1 = i.split('\t')[6].split('~')[1].split(':')[0]  # Get contig name1
            chm2 = i.split('\t')[6].split('~')[2].split(':')[0]  # Get contig name2
            col4 = i.split('\t')[3].split(' ')[1].split('~')[1]  # Get SV type label
            le = int(i.split('\t')[6].split('~')[1].split(':')[1])  # Get default left coord
            r = int(i.split('\t')[6].split('~')[2].split(':')[1])  # Get default right coord
            # Assign coord to left or right group depending on strand
            if col4.split(',')[0] == '+':
                leftclust[chm1][rnameidx + '-l'] = le
            elif col4.split(',')[0] == '-':
                rightclust[chm1][rnameidx + '-l'] = le
            if col4.split(',')[1] == '+':
                rightclust[chm2][rnameidx + '-r'] = r
            elif col4.split(',')[1] == '-':
                leftclust[chm2][rnameidx + '-r'] = r
        else:
            raise Exception('Error: %s is missed' % rnameidx)

    leftchrnamelist, leftchrcoorlist = OrderedDict(), OrderedDict()
    rightchrnamelist, rightchrcoorlist = OrderedDict(), OrderedDict()

    # Place the coords and read ids in a linear list for clustering with the same coord ascending order
    for chrm in contigs:
        leftclust[chrm] = OrderedDict(sorted(leftclust[chrm].items(), key=lambda y: y[1]))
        rightclust[chrm] = OrderedDict(sorted(rightclust[chrm].items(), key=lambda y: y[1]))
        leftchrnamelist[chrm] = []
        leftchrcoorlist[chrm] = []
        rightchrnamelist[chrm] = []
        rightchrcoorlist[chrm] = []
        for key in leftclust[chrm]:
            leftchrnamelist[chrm].append(key[:-2])
            leftchrcoorlist[chrm].append(leftclust[chrm][key])
        for key in rightclust[chrm]:
            rightchrnamelist[chrm].append(key[:-2])
            rightchrcoorlist[chrm].append(rightclust[chrm][key])

    # Carry out clustering
    readteam, mainclass = cluster(leftchrnamelist, leftchrcoorlist, rightchrnamelist, rightchrcoorlist, buf, svsizedict,
                                  classdict, ins_switch, mincov)
    """
    - readteam: Dictionary where keys are cluster coordinates (e.g. chr1:100-chr1-200 or chr1:100 for Ins SVs) and values are 
                lists of reads with the respective cluster coordinates, with the first read in the list as the lead read.
    - infodict: Dictionary of default data for each SV read from total_out
    - classdict: Dictionary of the SV class of each read (key=read_id)
    - mainclass: Dictionary where keys are cluster coordinates and values are the SV class of the lead read
    - svsizedict: Dictionary of SV sizes for each SV (key=read_id)
    """

    return readteam, infodict, classdict, mainclass, svsizedict


# Cluster coordinates
def cluster(leftchrnamelist,
            leftchrcoorlist,
            rightchrnamelist,
            rightchrcoorlist,
            buf,
            svsizedict,
            classdict,
            ins_switch,
            mincov):

    leftcluster2read, rightcluster2read, read2cluster = OrderedDict(), OrderedDict(), OrderedDict()

    # Cluster coordinates for each contig using buf as reference, clustering separately for left and right groups
    # Note: Coordinates within buf bp of each other will be clustered together and averaged to obtain a final coord. This is
    # done regardless of SV type or validity. Hence, final clustered coord may be slightly different from coord given by each
    # read.
    # Left group clustering
    for chrm in leftchrnamelist:
        n = len(leftchrnamelist[chrm])
        last = None
        group = []
        groupnames = []
        for i in range(n):
            if last is None or leftchrcoorlist[chrm][i] - last <= buf:
                group.append(leftchrcoorlist[chrm][i])
                groupnames.append(leftchrnamelist[chrm][i])
            else:
                avgcoord = int(round(sum(group) / len(group), 0))
                leftcluster2read[chrm + ':' + str(avgcoord)] = sorted(set(groupnames))
                for read in sorted(set(groupnames)):
                    if read not in read2cluster:
                        read2cluster[read] = []
                    read2cluster[read].append(chrm + ':' + str(avgcoord))
                group = [leftchrcoorlist[chrm][i]]
                groupnames = [leftchrnamelist[chrm][i]]
            last = leftchrcoorlist[chrm][i]
        if group:
            avgcoord = int(round(sum(group) / len(group), 0))
            leftcluster2read[chrm + ':' + str(avgcoord)] = sorted(set(groupnames))
            for read in sorted(set(groupnames)):
                if read not in read2cluster:
                    read2cluster[read] = []
                read2cluster[read].append(chrm + ':' + str(avgcoord))

    # Right group clustering
    for chrm in rightchrnamelist:
        n = len(rightchrnamelist[chrm])
        last = None
        group = []
        groupnames = []
        for i in range(n):
            if last is None or rightchrcoorlist[chrm][i] - last <= buf:
                group.append(rightchrcoorlist[chrm][i])
                groupnames.append(rightchrnamelist[chrm][i])
            else:
                avgcoord = int(round(sum(group) / len(group), 0))
                rightcluster2read[chrm + ':' + str(avgcoord)] = sorted(set(groupnames))
                for read in sorted(set(groupnames)):
                    if read not in read2cluster:
                        read2cluster[read] = []
                    read2cluster[read].append(chrm + ':' + str(avgcoord))
                group = [rightchrcoorlist[chrm][i]]
                groupnames = [rightchrnamelist[chrm][i]]
            last = rightchrcoorlist[chrm][i]
        if group:
            avgcoord = int(round(sum(group) / len(group), 0))
            rightcluster2read[chrm + ':' + str(avgcoord)] = sorted(set(groupnames))
            for read in sorted(set(groupnames)):
                if read not in read2cluster:
                    read2cluster[read] = []
                read2cluster[read].append(chrm + ':' + str(avgcoord))
    """
    - leftcluster2read or rightcluster2read: Dict, key=cluster, value=list of read_ids having the cluster
    - read2cluster: Dict, key=read_id, value=cluster
    """

    # Links each of the two coord in an SV together by read information (Only if SV has two coord, else coord will link to [])
    clusterlink = OrderedDict()
    coordlist2 = []

    for read in read2cluster:
        if read2cluster[read][0] not in clusterlink:  # Might replace with defaultdict in the future
            clusterlink[read2cluster[read][0]] = []
        if len(read2cluster[read]) == 1:  # If read only has single coord e.g. S-Nov_Ins or E-Nov_Ins, pass.
            pass
        elif len(read2cluster[read]) == 2:  # Append second coord to list
            clusterlink[read2cluster[read][0]].append(read2cluster[read][1])
            coordlist2.append(read2cluster[read][1])
        else:
            raise Exception("Read %s has none or more than 2 clusters" % read)

    # Filter 1
    # Clusters which have been regarded as second clusters will be removed as first clusters in dict to avoid double entry.
    for clusters in coordlist2:
        try:  # If second cluster is also a first cluster in clusterlink
            if not clusterlink[clusters]:  # If second cluster as a first cluster does not have own second cluster (single coord)
                del clusterlink[clusters]  # delete the cluster
        except KeyError:
            pass

    # Order second cluster lists chronologically
    for clusters in clusterlink:
        clusterlink[clusters] = sorted(set(clusterlink[clusters]))

    readteam, mainclass = OrderedDict(), OrderedDict()
    pastbps = {}

    # Assign the lead read and gather read team of each paired or unpaired cluster
    for coord in clusterlink:  # For clusters that are either unpaired or paired (left)
        if not clusterlink[coord]:  # If cluster is unpaired, which means single breakend SVs
            try:  # If coord in left group
                lead, main = leadread_bp(leftcluster2read[coord], svsizedict, classdict)  # Gather lead read
                readteam[coord] = [lead]
                readteam[coord].extend(sorted(set(leftcluster2read[coord]).difference([lead])))  # Append remaining reads
            except KeyError:  # If coord in right group
                lead, main = leadread_bp(rightcluster2read[coord], svsizedict, classdict)  # Gather lead read
                readteam[coord] = [lead]
                readteam[coord].extend(sorted(set(rightcluster2read[coord]).difference([lead])))  # Append remaining reads
            mainclass[coord] = main
        else:  # If cluster is paired (left)
            for coord2 in clusterlink[coord]:  # For each coord2 linked to coord1
                try:  # If coord1 is a cluster in left group
                    reads1 = leftcluster2read[coord]
                    try:  # If coord2 is a cluster in right group
                        reads2 = rightcluster2read[coord2]
                        if len(set(reads1).intersection(reads2)) < 1:  # Make sure the linked clusters are truly linked by >1 read
                            raise ValueError
                    except (KeyError, ValueError):  # If coord2 is a cluster in left group
                        reads2 = leftcluster2read[coord2]
                        if len(set(reads1).intersection(reads2)) < 1:  # Make sure the linked clusters are truly linked by >1 read
                            raise ValueError
                except (KeyError, ValueError):  # If coord1 is a cluster in right group
                    reads1 = rightcluster2read[coord]
                    try:   # If coord2 is a cluster in right group
                        reads2 = rightcluster2read[coord2]
                        if len(set(reads1).intersection(reads2)) < 1:  # Make sure the linked clusters are truly linked by >1 read
                            raise ValueError
                    except (KeyError, ValueError):  # If coord2 is a cluster in left group
                        reads2 = leftcluster2read[coord2]
                        if len(set(reads1).intersection(reads2)) < 1:  # Make sure the linked clusters are truly linked by >1 read
                            # No correct permutations for correct linked clusters
                            raise Exception('Error: Please raise this issue on GitHub, Thanks.')

                readlist = sorted(set(reads1).intersection(reads2))  # Get reads that intersect with both clusters
                for read in sorted(set(reads1).difference(reads2)) + sorted(set(reads2).difference(reads1)):
                    # If read only found in either clusters, include them only if they are bp_Nov_Ins
                    if classdict[read] == 'bp_Nov_Ins':
                        if read not in pastbps:
                            readlist.append(read)
                            pastbps[read] = 1

                # Gather lead read
                lead, main = leadread(sorted(set(readlist)), svsizedict, classdict, ins_switch)
                # For Nov_Ins, get the average coord between the two
                if main == 'Nov_Ins':
                    avg_coord = coord.split(':')[0] + ':' + \
                                str(int(round((int(coord.split(':')[1]) + int(coord2.split(':')[1])) / 2, 0)))
                    readteam[avg_coord] = [lead]
                    readteam[avg_coord].extend(sorted(set(readlist).difference([lead])))    # Append remaining reads
                    mainclass[avg_coord] = main
                else:  # Else include both coord into cluster label
                    readteam[coord + '-' + coord2] = [lead]
                    readteam[coord + '-' + coord2].extend(sorted(set(readlist).difference([lead])))    # Append remaining reads
                    mainclass[coord + '-' + coord2] = main

    # Filter 2
    # Delete cluster if it is supported by <mincov number of reads
    to_del = []
    for key, reads in readteam.items():
        if len(reads) < mincov:
            to_del.append(key)
    for key in to_del:
        del readteam[key]

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


def leadread_bp(reads, svsizedict, classdict):
    mainsvclass = 'bp_Nov_Ins'
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
        raise KeyError("Error: Main SV class not found")
    return leader, mainsvclass


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
        a = ','.join([y for y in tmpdict])
        raise Exception("Error: Unknown SV class %s, %s" % (a, ','.join(reads)))
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


# Function to generate BED from subdata
def normalbed(subdata):
    totalbed = []
    len_buf = 100  # Arbituary length deduction buffer for normal reads,
    for line in subdata:
        if int(line.split('\t')[1]) > 200:  # only allow len_buf if read start is 200 bp away from contig start
            coord1 = int(line.split('\t')[1]) + len_buf
        else:
            coord1 = int(line.split('\t')[1])
        # note that coord2 might be close to contig end but cannot be adjusted as in coord1 due to unknown contig length.
        coord2 = int(line.split('\t')[1]) + int(line.split('\t')[2]) - len_buf
        if coord2 - coord1 > 0:
            totalbed.append(line.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + line.split('\t')[4])
        else:
            pass
    return totalbed


# Create BED from SV breakends
def svbed2(readteam, mainclass):
    totalsvbed = []
    for clusters in readteam:
        svtype = mainclass[clusters]
        c = clusters.split('-')
        totalsvbed.append(c[0].split(':')[0] + '\t' + c[0].split(':')[1] + '\t' + str(int(c[0].split(':')[1]) + 1) + '\t' +
                          clusters + '\t' + svtype + '\t' + '1')
        if svtype in ['Nov_Ins', 'bp_Nov_Ins']:
            pass
        else:
            totalsvbed.append(c[1].split(':')[0] + '\t' + c[1].split(':')[1] + '\t' + str(int(c[1].split(':')[1]) + 1) + '\t' +
                              clusters + '\t' + svtype + '\t' + '2')
    return totalsvbed


# Parse BED intersection output
def intersection3(bed, readteam, newdictnouniq, mainclass):
    norm_connect = defaultdict(set)
    norm_cov = defaultdict(int)
    for line in bed:
        clusters = line.fields[7]
        ind = line.fields[9]
        read = line.fields[3]
        if read not in newdictnouniq[clusters]:
            norm_connect[clusters + '~' + ind].add(read)
    for clusters in readteam:
        svtype = mainclass[clusters]
        if svtype in ['Del', 'Inv', 'Inv2']:
            norm_cov[clusters] = int(round((len(norm_connect[clusters + '~1']) + len(norm_connect[clusters + '~2']))/2, 0))
        elif svtype == 'TDupl':
            norm_cov[clusters] = len(norm_connect[clusters + '~1'].intersection(norm_connect[clusters + '~2']))
        elif svtype in ['Nov_Ins', 'bp_Nov_Ins']:
            norm_cov[clusters] = len(norm_connect[clusters + '~1'])
        elif svtype in ['Intra-Ins2', 'Intra-Ins', 'Inter', 'Inter-Ins']:
            norm_cov[clusters] = min(len(norm_connect[clusters + '~1']), len(norm_connect[clusters + '~2']))
        else:
            raise Exception('Error: SV type not recognised.')
    return norm_cov


# Remove unique identifier from read name
def remove_uniq(readteam):
    newdictnouniq = {}
    for clusters in readteam:
        newdictnouniq[clusters] = set()
        for read in readteam[clusters]:
            newdictnouniq[clusters].add(read[:-6])
    return newdictnouniq


# Parse info into output
def arrange(svnormalcov, readteam, maxovl, mincov, infodict, mainclass, svsizedict):
    output = []
    n = 1
    for clusters in readteam:
        svtype = mainclass[clusters]
        lcov = countcov(readteam[clusters])
        bestread = readteam[clusters][0]
        # Filter 3
        if mincov <= lcov <= maxovl and svnormalcov[clusters] <= maxovl:
            if svtype in ['Inter', 'Inter-Ins']:
                _cluster = clusterorder(clusters,
                                        infodict[bestread].split('\t')[6].split('~')[1].split(':')[0],
                                        infodict[bestread].split('\t')[6].split('~')[1].split(':')[1],
                                        'Inter')
                output.append(
                    '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                    infodict[bestread].split('\t')[6].split('~')[0] +
                    '~' + _cluster.split('-')[0] + '~' + _cluster.split('-')[1] + '\t' +
                    '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                    ','.join(dotter(readteam[clusters])) + '\t' + str(svnormalcov[clusters])
                )
            elif svtype in ['Intra-Ins2', 'Intra-Ins']:
                _cluster = clusterorder(clusters,
                                        infodict[bestread].split('\t')[6].split('~')[1].split(':')[0],
                                        infodict[bestread].split('\t')[6].split('~')[1].split(':')[1].split('-')[0],
                                        'Intra')
                chrm_coord1 = _cluster.split('-')[0]
                coord2 = _cluster.split('-')[1].split(':')[1]
                output.append(
                    '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                    infodict[bestread].split('\t')[6].split('~')[0] + '~' +
                    chrm_coord1 + '-' + coord2 + '\t' +
                    '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                    ','.join(dotter(readteam[clusters])) + '\t' + str(svnormalcov[clusters])
                    )
            # Filter 3
            elif svtype in ['Nov_Ins', 'Del'] and svsizedict[bestread] <= 200 and lcov < max(mincov, 2):
                continue
            # Filter 4
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
                    chrm = clusters.split('-')[0].split(':')[0]
                    coord1 = str(min(int(clusters.split('-')[0].split(':')[1]), int(clusters.split('-')[1].split(':')[1])))
                    coord2 = str(max(int(clusters.split('-')[0].split(':')[1]), int(clusters.split('-')[1].split(':')[1])))
                    output.append(
                        '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                        infodict[bestread].split('\t')[6].split('~')[0] + '~' +
                        chrm + ':' + coord1 + '-' + coord2 + '\t' +
                        '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                        ','.join(dotter(readteam[clusters])) + '\t' + str(svnormalcov[clusters])
                    )
                else:
                    raise Exception('Error: Cluster name error')
            n += 1
    return output


# Function to convert non-digits to zero
def alpha(n):
    if str(n).isdigit():
        return int(n)
    else:
        return int(0)


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


# Convert single element list to .
def dotter(readlist):
    if len(readlist) == 1:
        return '.'
    else:
        return readlist[1:]


# Function to count total number of unique reads in an overlaping breakpoint entry
def countcov(readlist):
    newlist = []
    for read in readlist:
        newlist.append(read[:-6])
    return len(set(newlist))


# Order coord in cluster label according to default order
def clusterorder(_cluster, d_chrm1, d_coord1, mode):
    chrm1 = _cluster.split('-')[0].split(':')[0]
    coord1 = _cluster.split('-')[0].split(':')[1]
    chrm2 = _cluster.split('-')[1].split(':')[0]
    coord2 = _cluster.split('-')[1].split(':')[1]
    if mode == 'Intra':
        coords = [coord1, coord2]  # make cluster coords in to list
        index1 = (min([0, 1], key=lambda i: abs(int(coords[i]) - int(d_coord1))))
        if index1 == 0:
            index2 = 1
        else:
            index2 = 0
        new_cluster = chrm1 + ':' + str(coords[index1]) + '-' + chrm2 + ':' + str(coords[index2])
    elif mode == 'Inter':
        if chrm1 == chrm2:  # if same chrm
            coords = [coord1, coord2]
            index1 = (min(range(len(coords)), key=lambda i: abs(int(coords[i]) - int(d_coord1))))
            if index1 == 0:
                index2 = 1
            else:
                index2 = 0
            new_cluster = chrm1 + ':' + str(coords[index1]) + '-' + chrm2 + ':' + str(coords[index2])
        else:  # if chrm are different
            if chrm1 == d_chrm1:
                new_cluster = _cluster
            else:  # if chrm2 == d_chrm1
                new_cluster = _cluster.split('-')[1] + '-' + _cluster.split('-')[0]
    else:
        raise Exception('Error: wrong mode')
    return new_cluster


# Note: Some bps seen in parse_file would be missing
# in cluster_file because they are overlapped by another
# bp within the same read

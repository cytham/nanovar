"""
Functions for SV breakend clustering.

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
from collections import OrderedDict, defaultdict
from pybedtools import BedTool
from statistics import median


# SV clustering main function
def sv_cluster(
        subdata,
        parse,
        buf,
        maxovl,
        mincov,
        contigs,
        hsb_switch,
        seed
):
    """Returns a dictionary of clusters to SV breakends (mm) or list of clustered SV breakends in strings (hsb)

Keyword arguments:
- subdata: BED information for all reads in input
- parse: parse data from breakpoint_parser
- buf: coordinate clustering buffer in base pairs
- maxovl: maximum number of read overlap allowed for breakend-supporting or opposing- reads
- mincov: minimum number of breakend-supporting reads required for a breakend
- contigs: list of contig names from input reference genome
- hsb_switch: boolean to switch between pre and post blast re-evaluation
- seed: seed for SV index numbering

    """
    # Parse SV locations and cluster them
    # If not blast alignment
    if not hsb_switch:
        cluster_parse = rangecollect(parse,
                                     buf,
                                     contigs,
                                     hsb_switch,
                                     mincov
                                     )
        return cluster_parse
    else:  # If blast alignment
        clustid2coord, clustid2reads, clustid2class, clustid2ins_size, infodict, classdict, svsizedict = rangecollect(parse,
                                                                                                                      buf,
                                                                                                                      contigs,
                                                                                                                      hsb_switch,
                                                                                                                      mincov
                                                                                                                      )

        # Convert subdata to BED format
        bed3 = normalbed(subdata)
        bed3 = BedTool('\n'.join(bed3), from_string=True)
        bed3 = bed3.sort()

        # Convert SV breakends to BED format
        # bed4 = svbed2(readteam, mainclass)
        bed4 = svbed2(clustid2reads, clustid2class, clustid2coord)
        bed4 = BedTool('\n'.join(bed4), from_string=True)
        bed4 = bed4.sort()

        # Intersect all read regions to SV breakends
        intersect = bed3.intersect(bed4, wa=True, wb=True)

        # Create a duplicated readteam dict where reads have their unique identifier removed
        # newdictnouniq = remove_uniq(readteam)
        newdictnouniq = remove_uniq(clustid2reads)

        # Parse intersection BED to discover breakend-opposing reads for each breakend
        svnormalcov = intersection3(intersect, clustid2reads, newdictnouniq, clustid2class)

        # Output all information into a list of strings
        output, seed = arrange(svnormalcov, clustid2coord, clustid2reads, maxovl, mincov, infodict, clustid2class, svsizedict, seed, clustid2ins_size)
        return output, seed + 1


# Function to collect information of each breakpoint entry and cluster breakpoints
def rangecollect(x, buf, contigs, hsb_switch, mincov):

    classdict, svsizedict, infodict = OrderedDict(), OrderedDict(), OrderedDict()
    leftclust, rightclust = OrderedDict(), OrderedDict()

    for c in contigs:
        leftclust[c] = OrderedDict()
        rightclust[c] = OrderedDict()

    for i in x:
        rnameidx = i.split('\t')[8]  # Read unique id
        infodict[rnameidx] = '\t'.join(i.split('\t')[0:10])  # Read info storage from total_out file
        svtype = str(svtypecorrector(i.split('\t')[3].split(' ')[0]))  # Get sv type and slightly modify it
        classdict[rnameidx] = svtype  # Assign sv type to classdict
        svsize = alpha(i.split('\t')[3].split(' ')[1].split('~')[0])  # Get sv size
        svsizedict[rnameidx] = svsize  # Assign sv size to svsizedict
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
                    raise Exception('Error: single SV breakend %s has unknown naming' % rnameidx)
            elif len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 2:  # If double breakends
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]  # Get contig name
                le = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])  # Get default left coord
                r = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[1])  # Get default right coord
                if svtype == 'TDupl':  # If tandem dup, reverse left and right groups
                    leftclust[chm1][rnameidx + '-r'] = r
                    rightclust[chm1][rnameidx + '-l'] = le
                elif svtype in ['Inv', 'Inv2']:
                    col4s = i.split('\t')[3].split(' ')[1].split('~')[1]  # Get SV strand info
                    if col4s.split(',')[0] == '+' and col4s.split(',')[1] == '-':
                        leftclust[chm1][rnameidx + '-l'] = le
                        leftclust[chm1][rnameidx + '-r'] = r
                    elif col4s.split(',')[0] == '-' and col4s.split(',')[1] == '+':
                        rightclust[chm1][rnameidx + '-l'] = le
                        rightclust[chm1][rnameidx + '-r'] = r
                    elif col4s.split(',')[0] == '+' and col4s.split(',')[1] == '+':
                        leftclust[chm1][rnameidx + '-l'] = le
                        rightclust[chm1][rnameidx + '-r'] = r
                    elif col4s.split(',')[0] == '-' and col4s.split(',')[1] == '-':
                        rightclust[chm1][rnameidx + '-l'] = le
                        leftclust[chm1][rnameidx + '-r'] = r
                    else:
                        raise Exception('Error: Inversion % has unidentified strand orientation' % rnameidx)
                elif svtype in ['Intra-Ins', 'Intra-Ins2']:
                    col4s = i.split('\t')[3].split(' ')[1].split('~')[1]  # Get SV strand info
                    if col4s.split(',')[0] == '+':
                        leftclust[chm1][rnameidx + '-l'] = le
                    elif col4s.split(',')[0] == '-':
                        rightclust[chm1][rnameidx + '-l'] = le
                    if col4s.split(',')[1] == '+':
                        rightclust[chm1][rnameidx + '-r'] = r
                    elif col4s.split(',')[1] == '-':
                        leftclust[chm1][rnameidx + '-r'] = r
                else:
                    leftclust[chm1][rnameidx + '-l'] = le
                    rightclust[chm1][rnameidx + '-r'] = r
        elif len(i.split('\t')[6].split('~')) == 3:  # Inter translocation
            chm1 = i.split('\t')[6].split('~')[1].split(':')[0]  # Get contig name1
            chm2 = i.split('\t')[6].split('~')[2].split(':')[0]  # Get contig name2
            col4s = i.split('\t')[3].split(' ')[1].split('~')[1]  # Get SV strand info
            le = int(i.split('\t')[6].split('~')[1].split(':')[1])  # Get default left coord
            r = int(i.split('\t')[6].split('~')[2].split(':')[1])  # Get default right coord
            # Assign coord to left or right group depending on strand
            if col4s.split(',')[0] == '+':
                leftclust[chm1][rnameidx + '-l'] = le
            elif col4s.split(',')[0] == '-':
                rightclust[chm1][rnameidx + '-l'] = le
            if col4s.split(',')[1] == '+':
                rightclust[chm2][rnameidx + '-r'] = r
            elif col4s.split(',')[1] == '-':
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
            leftchrnamelist[chrm].append(key)  # [:-2]
            leftchrcoorlist[chrm].append(leftclust[chrm][key])
        for key in rightclust[chrm]:
            rightchrnamelist[chrm].append(key)  # [:-2]
            rightchrcoorlist[chrm].append(rightclust[chrm][key])
    
    # Carry out clustering
    if not hsb_switch:
        cluster_parse = cluster(leftchrnamelist, leftchrcoorlist, rightchrnamelist, rightchrcoorlist, buf, svsizedict,
                                classdict, hsb_switch, mincov)
        return cluster_parse
    else:
        clustid2coord, clustid2reads, clustid2class, clustid2ins_size = cluster(leftchrnamelist, leftchrcoorlist, rightchrnamelist, rightchrcoorlist, buf, svsizedict,
                                                                        classdict, hsb_switch, mincov)
        return clustid2coord, clustid2reads, clustid2class, clustid2ins_size, infodict, classdict, svsizedict


"""
- cluster_parse: Dictionary where keys are the cluster coordinates + svtype and values are the breakend names respective to each 
                cluster.
- readteam: Dictionary where keys are cluster coordinates (e.g. chr1:100-chr1-200 or chr1:100 for Ins SVs) and values are 
            lists of reads with the respective cluster coordinates, with the first read in the list as the lead read.
- infodict: Dictionary of default data for each SV read from total_out
- classdict: Dictionary of the SV class of each read (key=read_id)
- mainclass: Dictionary where keys are cluster coordinates and values are the SV class of the lead read
- svsizedict: Dictionary of SV sizes for each SV (key=read_id)
"""


# Cluster coordinates
def cluster(leftchrnamelist,
            leftchrcoorlist,
            rightchrnamelist,
            rightchrcoorlist,
            buf,
            svsizedict,
            classdict,
            hsb_switch,
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
                # avgcoord = int(round(sum(group) / len(group), 0))
                medcoord = int(median(group))
                leftcluster2read[chrm + ':' + str(medcoord)] = sorted(set(y[:-2] for y in groupnames))
                for read in sorted(set(groupnames)):
                    if read[:-2] not in read2cluster:
                        read2cluster[read[:-2]] = []
                    read2cluster[read[:-2]].append(chrm + ':' + str(medcoord) + 'l')
                group = [leftchrcoorlist[chrm][i]]
                groupnames = [leftchrnamelist[chrm][i]]
            last = leftchrcoorlist[chrm][i]
        if group:
            # avgcoord = int(round(sum(group) / len(group), 0))
            medcoord = int(median(group))
            leftcluster2read[chrm + ':' + str(medcoord)] = sorted(set(y[:-2] for y in groupnames))
            for read in sorted(set(groupnames)):
                if read[:-2] not in read2cluster:
                    read2cluster[read[:-2]] = []
                read2cluster[read[:-2]].append(chrm + ':' + str(medcoord) + 'l')
    
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
                # avgcoord = int(round(sum(group) / len(group), 0))
                medcoord = int(median(group))
                rightcluster2read[chrm + ':' + str(medcoord)] = sorted(set(y[:-2] for y in groupnames))
                for read in sorted(set(groupnames)):
                    if read[:-2] not in read2cluster:
                        read2cluster[read[:-2]] = []
                    read2cluster[read[:-2]].append(chrm + ':' + str(medcoord) + 'r')
                group = [rightchrcoorlist[chrm][i]]
                groupnames = [rightchrnamelist[chrm][i]]
            last = rightchrcoorlist[chrm][i]
        if group:
            # avgcoord = int(round(sum(group) / len(group), 0))
            medcoord = int(median(group))
            rightcluster2read[chrm + ':' + str(medcoord)] = sorted(set(y[:-2] for y in groupnames))
            for read in sorted(set(groupnames)):
                if read[:-2] not in read2cluster:
                    read2cluster[read[:-2]] = []
                read2cluster[read[:-2]].append(chrm + ':' + str(medcoord) + 'r')
    """
    - leftcluster2read or rightcluster2read: Dict, key=cluster, value=list of read_ids having the cluster
    - read2cluster: Dict, key=read_id, value=cluster
    """
    
    # Links each of the two coord in an SV together by read information (Only if SV has two coord, else coord will link to [])
    clusterlink = OrderedDict()
    coordlist2 = []
    for read in read2cluster:
        if read2cluster[read][0] not in clusterlink:
            clusterlink[read2cluster[read][0]] = []
        if len(read2cluster[read]) == 1:  # If read only has single coord e.g. S-Nov_Ins or E-Nov_Ins, pass.
            pass
        elif len(read2cluster[read]) == 2:  # Append second coord to list
            clusterlink[read2cluster[read][0]].append(read2cluster[read][1])
            coordlist2.append(read2cluster[read][1])
        else:
            raise Exception("Error: Read %s has %d clusters" % (read, len(read2cluster[read])))
            # Read has 0 or more than 2 clusters

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
    
    # readteam, mainclass = OrderedDict(), OrderedDict()
    # pastbps = {}
    # cluster_parse = {}
    # ins_median_size = {}
    # for coord in clusterlink:  # For clusters that are either unpaired or paired (left)
    #     if not clusterlink[coord]:  # If cluster is unpaired, which means single breakend SVs
    #         if not hsb_switch:
    #             # if coord[-1] == 'l':
    #             #     readlist = leftcluster2read[coord[:-1]]
    #             # elif coord[-1] == 'r':
    #             #     readlist = rightcluster2read[coord[:-1]]
    #             # else:
    #             #     raise Exception('Error: Coord %s did not specify direction' % coord)
    #             # lead, main = leadread_bp(readlist, svsizedict, classdict)  # Gather lead read
    #             # if len({x.rsplit('~', 1)[0] for x in readlist}) >= min(2, mincov):
    #             #     cluster_parse[coord + '~' + main] = sorted(readlist)
    #             pass
    #         else:
    #             if coord[-1] == 'l':
    #                 readlist = leftcluster2read[coord[:-1]]
    #             elif coord[-1] == 'r':
    #                 readlist = rightcluster2read[coord[:-1]]
    #             else:
    #                 raise Exception('Error: Coord %s did not specify direction' % coord)
    #             # mm = []
    #             # hsb = []
    #             # for r in readlist:
    #             #     if len(r.rsplit('~', 1)[1]) == 5:  # mm
    #             #         mm.append(r)
    #             #     elif len(r.rsplit('~', 1)[1]) == 4:  # hsb
    #             #         hsb.append(r.rsplit('~', 1)[0])
    #             #     else:
    #             #         raise Exception('Error: %s read has unrecognised unique ID' % r)
    #             # for r in mm:
    #             #     if r.rsplit('~', 1)[0] in hsb:
    #             #         readlist.remove(r)  # Remove minimap reads if blast reads present in read list
    #             lead, main, med_size = leadread_bp(readlist, svsizedict, classdict)  # Gather lead read
    #             readteam[coord] = [lead]
    #             readteam[coord].extend(sorted(set(readlist).difference([lead])))  # Append remaining reads
    #             mainclass[coord] = main
    #             ins_median_size[coord] = med_size
    #     else:  # If cluster is paired (left)
    #         if coord[-1] == 'l':
    #             reads1 = leftcluster2read[coord[:-1]]
    #         elif coord[-1] == 'r':
    #             reads1 = rightcluster2read[coord[:-1]]
    #         else:
    #             raise Exception('Error: Coord %s did not specify direction' % coord)
    #         for coord2 in clusterlink[coord]:
    #             if coord2[-1] == 'l':
    #                 reads2 = leftcluster2read[coord2[:-1]]
    #             elif coord2[-1] == 'r':
    #                 reads2 = rightcluster2read[coord2[:-1]]
    #             else:
    #                 raise Exception('Error: Coord %s did not specify direction' % coord2)
    #             if coord[-1] in ['l', 'r']:
    #                 coord = coord[:-1]
    #             coord2 = coord2[:-1]
    #             #readlist0 = sorted(set(reads1).intersection(reads2))  # Get reads that intersect with both clusters
    #             readlists = separate_cluster_type(reads1, reads2, classdict)                
    #             if readlists:
    #                 for _readlist in readlists:
    #                     readlist = [x for x in _readlist if x not in pastbps]
    #                     if readlist:
    #                         for i in readlist:
    #                             pastbps[i] = 1
    #                         if not hsb_switch:
    #                             # for read in sorted(set(reads1).difference(reads2)) + sorted(set(reads2).difference(reads1)):
    #                             #     if classdict[read] == 'bp_Nov_Ins':
    #                             #         if read not in pastbps:
    #                             #             readlist.append(read)
    #                             #             pastbps[read] = 1
    #                             # lead, main = leadread(sorted(set(readlist)), svsizedict, classdict, hsb_switch)
    #                             # # Filter 2
    #                             # # Delete cluster if it is supported by <mincov or 2 number of reads, else if it is <100 bp
    #                             # if main == 'Del':
    #                             #     if len({x.rsplit('~', 1)[0] for x in readlist}) >= min(2, mincov):
    #                             #         cluster_parse[coord + '-' + coord2 + '~' + main] = sorted(readlist)
    #                             #     else:
    #                             #         if svsizedict[readlist[0]] >= 100:  # Arbituary threshold for small deletions
    #                             #             cluster_parse[coord + '-' + coord2 + '~' + main] = sorted(readlist)
    #                             # elif main == 'Nov_Ins':
    #                             #     if len({x.rsplit('~', 1)[0] for x in readlist}) >= min(2, mincov):
    #                             #         avg_coord = coord.split(':')[0] + ':' + \
    #                             #                     str(int(round((int(coord.split(':')[1]) + int(coord2.split(':')[1])) / 2, 0)))
    #                             #         cluster_parse[avg_coord + '~' + main] = sorted(readlist)
    #                             #     else:
    #                             #         if svsizedict[readlist[0]] >= 100:  # Arbituary threshold for small insertions
    #                             #             avg_coord = coord.split(':')[0] + ':' + \
    #                             #                         str(int(round((int(coord.split(':')[1]) + int(coord2.split(':')[1])) / 2, 0)))
    #                             #             cluster_parse[avg_coord + '~' + main] = sorted(readlist)
    #                             # else:
    #                             #     cluster_parse[coord + '-' + coord2 + '~' + main] = sorted(readlist)
    #                             pass
    #                         else:
    #                             lead, main, med_size = leadread(sorted(set(readlist)), svsizedict, classdict, hsb_switch)
    #                             # Check for reads exclusive to each left/right cluster and add them into readlist if conditions satisfied
    #                             for read in sorted(set(reads1).difference(reads2)) + sorted(set(reads2).difference(reads1)):
    #                                 if classdict[read] == 'bp_Nov_Ins':
    #                                     if read not in pastbps:
    #                                         readlist.append(read)
    #                                         pastbps[read] = 1
    #                                 else:
    #                                     if main == 'Del':
    #                                         if classdict[read] in ['Del']:
    #                                             if read not in pastbps:
    #                                                 readlist.append(read)
    #                                                 pastbps[read] = 1
    #                                     elif main in ['Intra-Ins', 'Intra-Ins2']:
    #                                         if classdict[read] in ['Nov_Ins', 'Intra-Ins', 'Intra-Ins2']:
    #                                             if read not in pastbps:
    #                                                 readlist.append(read)
    #                                                 pastbps[read] = 1
    #                                     elif main in ['Inter', 'Inter-Ins']:
    #                                         if classdict[read] in ['Nov_Ins', 'Inter', 'Inter-Ins']:
    #                                             if read not in pastbps:
    #                                                 readlist.append(read)
    #                                                 pastbps[read] = 1
    #                                     elif main in ['Inv', 'Inv2']:
    #                                         if classdict[read] in ['Nov_Ins', 'Inv', 'Inv2']:
    #                                             if read not in pastbps:
    #                                                 readlist.append(read)
    #                                                 pastbps[read] = 1
    #                                     elif main in 'TDupl':
    #                                         if classdict[read] in ['Nov_Ins', 'TDupl']:
    #                                             if read not in pastbps:
    #                                                 readlist.append(read)
    #                                                 pastbps[read] = 1
    #                                     elif main == 'Nov_Ins':
    #                                         if classdict[read] in ['Intra-Ins', 'Intra-Ins2', 'Inter', 'Inter-Ins', 'TDupl']:
    #                                             if read not in pastbps:
    #                                                 readlist.append(read)
    #                                                 pastbps[read] = 1
    #                                     else:
    #                                         raise Exception('Error: Main SV type in cluster %s with lead %s is not recognised' %
    #                                                         ((coord+coord2), lead))
    #                             # Remove duplicate reads from mm
    #                             # mm = []
    #                             # hsb = []
    #                             # for r in readlist:
    #                             #     if len(r.rsplit('~', 1)[1]) == 5:  # mm
    #                             #         mm.append(r)
    #                             #     elif len(r.rsplit('~', 1)[1]) == 4:  # hsb
    #                             #         hsb.append(r.rsplit('~', 1)[0])
    #                             #     else:
    #                             #         raise Exception('Error: %d read has unrecognised unique ID' % r)
    #                             # for r in mm:
    #                             #     rn = r.rsplit('~', 1)[0]
    #                             #     if rn in hsb:
    #                             #         if r == lead:
    #                             #             # readlist.remove(r)  # Remove blast read if mm read is leader present in read list
    #                             #             for rl in readlist:
    #                             #                 if rl.rsplit('~', 1)[0] == rn and len(rl.rsplit('~', 1)[1]) == 4:
    #                             #                     readlist.remove(rl)
    #                             #                     # lead = rl
    #                             #                     # main = classdict[rl]
    #                             #         else:
    #                             #             readlist.remove(r)  # Remove minimap reads if blast reads present in read list
    #                             # lead, main = leadread(sorted(set(readlist)), svsizedict, classdict, hsb_switch)
    #                             if main == 'Nov_Ins':
    #                                 avg_coord = coord.split(':')[0] + ':' + \
    #                                             str(int(round((int(coord.split(':')[1]) + int(coord2.split(':')[1])) / 2, 0)))
    #                                 readteam[avg_coord] = [lead]
    #                                 readteam[avg_coord].extend(sorted(set(readlist).difference([lead])))  # Append remaining reads
    #                                 mainclass[avg_coord] = main
    #                                 ins_median_size[avg_coord] = med_size
    #                             else:  # Else include both coord into cluster label
    #                                 readteam[coord + '-' + coord2] = [lead]
    #                                 readteam[coord + '-' + coord2].extend(
    #                                     sorted(set(readlist).difference([lead])))  # Append remaining reads
    #                                 mainclass[coord + '-' + coord2] = main
    pastbps = {}
    clustid = 0
    clustid2coord = OrderedDict()
    clustid2reads = OrderedDict()
    clustid2class = OrderedDict()
    clustid2ins_size = OrderedDict()
    for coord in clusterlink:  # For clusters that are either unpaired or paired (left)
        if not clusterlink[coord]:  # If cluster is unpaired, which means single breakend SVs
            if coord[-1] == 'l':
                readlist = leftcluster2read[coord[:-1]]
            elif coord[-1] == 'r':
                readlist = rightcluster2read[coord[:-1]]
            else:
                raise Exception('Error: Coord %s did not specify direction' % coord)
            lead, main, med_size = leadread_bp(readlist, svsizedict, classdict)  # Gather lead read
            clustid2reads[clustid] = [lead]
            clustid2reads[clustid].extend(sorted(set(readlist).difference([lead])))  # Append remaining reads
            clustid2coord[clustid] = coord
            clustid2class[clustid] = main
            clustid2ins_size[clustid] = med_size
            clustid += 1
        else:  # If cluster is paired (left)
            if coord[-1] == 'l':
                reads1 = leftcluster2read[coord[:-1]]
            elif coord[-1] == 'r':
                reads1 = rightcluster2read[coord[:-1]]
            else:
                raise Exception('Error: Coord %s did not specify direction' % coord)
            for coord2 in clusterlink[coord]:
                if coord2[-1] == 'l':
                    reads2 = leftcluster2read[coord2[:-1]]
                elif coord2[-1] == 'r':
                    reads2 = rightcluster2read[coord2[:-1]]
                else:
                    raise Exception('Error: Coord %s did not specify direction' % coord2)
                if coord[-1] in ['l', 'r']:
                    coord = coord[:-1]
                coord2 = coord2[:-1]
                readlists, remainder = separate_cluster_type(reads1, reads2, classdict)                
                if readlists:
                    for _readlist in readlists:
                        readlist = [x for x in _readlist if x not in pastbps]
                        if readlist:
                            for i in readlist:
                                pastbps[i] = 1
                            lead, main, med_size = leadread(sorted(set(readlist)), svsizedict, classdict, hsb_switch)
                            # Check for reads exclusive to each left/right cluster and add them into readlist if conditions satisfied
                            for read in remainder:
                                if classdict[read] == 'bp_Nov_Ins':
                                    if read not in pastbps:
                                        readlist.append(read)
                                        pastbps[read] = 1
                                else:
                                    if main == 'Del':
                                        if classdict[read] in ['Del']:
                                            if read not in pastbps:
                                                readlist.append(read)
                                                pastbps[read] = 1
                                    elif main in ['Intra-Ins', 'Intra-Ins2']:
                                        if classdict[read] in ['Nov_Ins', 'Intra-Ins', 'Intra-Ins2']:
                                            if read not in pastbps:
                                                readlist.append(read)
                                                pastbps[read] = 1
                                    elif main in ['Inter', 'Inter-Ins']:
                                        if classdict[read] in ['Nov_Ins', 'Inter', 'Inter-Ins']:
                                            if read not in pastbps:
                                                readlist.append(read)
                                                pastbps[read] = 1
                                    elif main in ['Inv', 'Inv2']:
                                        if classdict[read] in ['Nov_Ins', 'Inv', 'Inv2']:
                                            if read not in pastbps:
                                                readlist.append(read)
                                                pastbps[read] = 1
                                    elif main in 'TDupl':
                                        if classdict[read] in ['Nov_Ins', 'TDupl']:
                                            if read not in pastbps:
                                                readlist.append(read)
                                                pastbps[read] = 1
                                    elif main == 'Nov_Ins':
                                        if classdict[read] in ['Intra-Ins', 'Intra-Ins2', 'Inter', 'Inter-Ins', 'TDupl']:
                                            if read not in pastbps:
                                                readlist.append(read)
                                                pastbps[read] = 1
                                    else:
                                        raise Exception('Error: Main SV type in cluster %s with lead %s is not recognised' %
                                                        ((coord+coord2), lead))
                            if main == 'Nov_Ins':
                                avg_coord = coord.split(':')[0] + ':' + \
                                            str(int(round((int(coord.split(':')[1]) + int(coord2.split(':')[1])) / 2, 0)))
                                clustid2reads[clustid] = [lead]
                                clustid2reads[clustid].extend(sorted(set(readlist).difference([lead])))  # Append remaining reads
                                clustid2coord[clustid] = avg_coord
                                clustid2class[clustid] = main
                                clustid2ins_size[clustid] = med_size
                                clustid += 1
                            else:  # Else include both coord into cluster label
                                clustid2reads[clustid] = [lead]
                                clustid2reads[clustid].extend(sorted(set(readlist).difference([lead])))  # Append remaining reads
                                clustid2coord[clustid] = coord + '-' + coord2
                                clustid2class[clustid] = main
                                clustid += 1
    # # Filter 2
    # # Delete cluster if it is supported by <mincov or 2 number of reads
    # if hsb_switch:
    #     to_del = []
    #     for key, reads in clustid2reads.items():
    #         if len({x.rsplit('~', 1)[0] for x in reads}) < min(2, mincov):
    #             to_del.append(key)
    #     # Delete all low cov clusters
    #     for key in to_del:
    #         del readteam[key]
    if not hsb_switch:
        #return cluster_parse
        pass
    else:
        return clustid2coord, clustid2reads, clustid2class, clustid2ins_size

# Separate clusters by per read svtype (this would deprecate mainclasssv)
def separate_cluster_type(reads1, reads2, classdict):
    intersect = sorted(set(reads1).intersection(reads2))
    cluster_class_dict = {}
    remainder = []
    for read in intersect:
        svtype = classdict[read][0:3]  # Taking first 3 element of svtype (E.g. Inv2 -> Inv)
        if svtype != 'bp_':
            if svtype not in cluster_class_dict:  
                cluster_class_dict[svtype] = [read]
            else:
                cluster_class_dict[svtype].append(read)
        else:
            remainder.append(read)
    _readlists = []
    for svtype in cluster_class_dict:
        if len(cluster_class_dict[svtype]) >=2:
            _readlists.append(cluster_class_dict[svtype])
    for read in sorted(set(reads1).difference(reads2)) + sorted(set(reads2).difference(reads1)):
        remainder.append(read)
    return _readlists, remainder

def leadread(reads, svsizedict, classdict, hsb_switch):
    mainsvclass = mainclasssv(reads, classdict, hsb_switch)
    sizedict = OrderedDict()
    leader = ''
    sv_sizes = []
    for read in reads:
        sizedict[read] = svsizedict[read]
        sv_sizes.append(svsizedict[read])
    sizedictsort = [key for (key, value) in sorted(sizedict.items(), key=lambda y: y[1], reverse=True)]
    for read in sizedictsort:
        if classdict[read] == mainsvclass:
            leader = read
            break
    if leader == '':
        raise Exception("Error: Main SV class not found")
    if mainsvclass == 'Nov_Ins':
        med_size = median(sv_sizes)
    else:
        med_size = 0
    return leader, mainsvclass, med_size


def leadread_bp(reads, svsizedict, classdict):
    mainsvclass = 'bp_Nov_Ins'
    sizedict = OrderedDict()
    leader = ''
    sv_sizes = []
    for read in reads:
        sizedict[read] = svsizedict[read]
        sv_sizes.append(svsizedict[read])
    sizedictsort = [key for (key, value) in sorted(sizedict.items(), key=lambda y: y[1], reverse=True)]
    for read in sizedictsort:
        if classdict[read] == mainsvclass:  # Read with the largest SV size and correspond to mainsvclass will be the leader
            leader = read
            break
    if leader == '':
        raise Exception("Error: Main SV class not found for reads %s" % ','.join(reads))
    med_size = median(sv_sizes)
    return leader, mainsvclass, med_size


# Function to rank and filter best sv type
def mainclasssv(reads, classdict, hsb_switch):
    if hsb_switch:
        tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
        tierxsv = ['Intra-Ins', 'Inter']
        tier2sv = ['Inv', 'Del', 'Nov_Ins', 'TDupl']
        tier1sv = ['bp_Nov_Ins']
    else:
        tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
        tierxsv = ['Intra-Ins', 'Inter']
        tier2sv = ['Inv', 'Del', 'TDupl', 'Nov_Ins']
        tier1sv = ['bp_Nov_Ins']
    tmpdict = defaultdict(int)
    for read in reads:
        svclass = classdict[read]
        tmpdict[svclass] += 1
    if any(y in tier3sv for y in tmpdict):
        for s in tierxsv:
            try:
                del tmpdict[s]
            except KeyError:
                pass
        for s in tier2sv:
            try:
                del tmpdict[s]
            except KeyError:
                pass
        try:
            del tmpdict["bp_Nov_Ins"]
        except KeyError:
            pass
    elif any(y in tierxsv for y in tmpdict):
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
def svbed2(clustid2reads, clustid2class, clustid2coord):
    totalsvbed = []
    for clustid in clustid2reads:
        svtype = clustid2class[clustid]
        c = [x.strip('lr') for x in clustid2coord[clustid].split('-')]
        totalsvbed.append(c[0].split(':')[0] + '\t' + c[0].split(':')[1] + '\t' + str(int(c[0].split(':')[1]) + 1) + '\t' +
                          str(clustid) + '\t' + svtype + '\t' + '1')
        if svtype in ['Nov_Ins', 'bp_Nov_Ins']:
            pass
        else:
            totalsvbed.append(c[1].split(':')[0] + '\t' + c[1].split(':')[1] + '\t' + str(int(c[1].split(':')[1]) + 1) + '\t' +
                              str(clustid) + '\t' + svtype + '\t' + '2')
    return totalsvbed


# Parse BED intersection output
def intersection3(bed, clustid2reads, newdictnouniq, clustid2class):
    norm_connect = defaultdict(set)
    norm_cov = defaultdict(int)
    for line in bed:
        clustid = line.fields[7]
        ind = line.fields[9]
        read = line.fields[3]
        if read not in newdictnouniq[int(clustid)]:
            norm_connect[clustid + '~' + ind].add(read)
    for clustid in clustid2reads:
        svtype = clustid2class[clustid]
        if svtype in ['Del', 'Inv', 'Inv2']:
            norm_cov[clustid] = int(round((len(norm_connect[str(clustid) + '~1']) + len(norm_connect[str(clustid) + '~2']))/2, 0))
        elif svtype == 'TDupl':
            norm_cov[clustid] = len(norm_connect[str(clustid) + '~1'].intersection(norm_connect[str(clustid) + '~2']))
        elif svtype in ['Nov_Ins', 'bp_Nov_Ins']:
            norm_cov[clustid] = len(norm_connect[str(clustid) + '~1'])
        elif svtype in ['Intra-Ins2', 'Intra-Ins', 'Inter', 'Inter-Ins']:
            norm_cov[clustid] = min(len(norm_connect[str(clustid) + '~1']), len(norm_connect[str(clustid) + '~2']))
        else:
            raise Exception('Error: SV type not recognised.')
    return norm_cov


# Remove unique identifier from read name
def remove_uniq(clustid2reads):
    newdictnouniq = {}
    for clustid in clustid2reads:
        newdictnouniq[clustid] = set()
        for read in clustid2reads[clustid]:
            newdictnouniq[clustid].add(read.rsplit('~', 1)[0])
    return newdictnouniq


# Parse info into output
def arrange(svnormalcov, clustid2coord, clustid2reads, maxovl, mincov, infodict, clustid2class, svsizedict, seed, clustid2ins_size):
    output = []
    n = seed
    for clustid in clustid2reads:
        svtype = clustid2class[clustid]
        lcov = countcov(clustid2reads[clustid])
        bestread = clustid2reads[clustid][0]
        # Filter 3
        if mincov <= lcov <= maxovl and svnormalcov[clustid] <= maxovl:
            if svtype in ['Inter', 'Inter-Ins']:
                _cluster = clusterorder(clustid2coord[clustid],
                                        infodict[bestread].split('\t')[6].split('~')[1].split(':')[0],
                                        infodict[bestread].split('\t')[6].split('~')[1].split(':')[1],
                                        'Inter')
                output.append(
                    '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                    infodict[bestread].split('\t')[6].split('~')[0] +
                    '~' + _cluster.split('-')[0] + '~' + _cluster.split('-')[1] + '\t' +
                    '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                    ','.join(dotter(clustid2reads[clustid])) + '\t' + str(svnormalcov[clustid])
                )
            elif svtype in ['Intra-Ins2', 'Intra-Ins']:
                _cluster = clusterorder(clustid2coord[clustid],
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
                    ','.join(dotter(clustid2reads[clustid])) + '\t' + str(svnormalcov[clustid])
                    )
            # Filter 3
            # elif svtype in ['Del'] and svsizedict[bestread] <= 200 and lcov < max(mincov, 2):  # minimally requires 2
            #     continue
            # elif svtype in ['Nov_Ins', 'bp_Nov_Ins'] and ins_median_size[clusters] <= 200 and lcov < max(mincov, 2):  # minimally requires 2
            #     continue
            else:
                if len(clustid2coord[clustid].split('-')) == 1:
                    output.append(
                        ins_size_mod(infodict[bestread].split('\t')[0:6], clustid2ins_size[clustid]) + '\tnv_SV' + str(n) + '-' +
                        infodict[bestread].split('\t')[6].split('~')[0] + '~' + clustid2coord[clustid].strip('lr') + '-' +
                        str(int(clustid2coord[clustid].strip('lr').split(':')[1]) + 1) + '\t' +
                        '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                        ','.join(dotter(clustid2reads[clustid])) + '\t' + str(svnormalcov[clustid])
                    )
                elif len(clustid2coord[clustid].split('-')) == 2:
                    chrm = clustid2coord[clustid].split('-')[0].split(':')[0]
                    coord1 = str(min(int(clustid2coord[clustid].split('-')[0].split(':')[1]), int(clustid2coord[clustid].split('-')[1].split(':')[1])))
                    coord2 = str(max(int(clustid2coord[clustid].split('-')[0].split(':')[1]), int(clustid2coord[clustid].split('-')[1].split(':')[1])))
                    if svtype == 'Nov_Ins':
                        output.append(
                            ins_size_mod(infodict[bestread].split('\t')[0:6], clustid2ins_size[clustid]) + '\tnv_SV' + str(n) + '-' +
                            infodict[bestread].split('\t')[6].split('~')[0] + '~' +
                            chrm + ':' + coord1 + '-' + coord2 + '\t' +
                            '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                            ','.join(dotter(clustid2reads[clustid])) + '\t' + str(svnormalcov[clustid])
                        )
                    else:
                        output.append(
                            '\t'.join(infodict[bestread].split('\t')[0:6]) + '\tnv_SV' + str(n) + '-' +
                            infodict[bestread].split('\t')[6].split('~')[0] + '~' +
                            chrm + ':' + coord1 + '-' + coord2 + '\t' +
                            '\t'.join(infodict[bestread].split('\t')[7:]) + '\t' + str(lcov) + '\t' +
                            ','.join(dotter(clustid2reads[clustid])) + '\t' + str(svnormalcov[clustid])
                        )
                else:
                    raise Exception('Error: Cluster name error')
            n += 1
    return output, n


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
        newlist.append(read.rsplit('~', 1)[0])
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


# Modify sv size in INS
def ins_size_mod(line, size):
    return '\t'.join(line[0:3] + [line[3].split(' ')[0] + ' ' + str(size) + '~' + line[3].split('~')[1]] + line[4:])
    

# Note: Some bps seen in parse_file would be missing
# in cluster_file because they are overlapped by another
# bp within the same read

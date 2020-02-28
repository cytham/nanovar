import sys
from sys import argv
from pybedtools import BedTool


if len(argv) != 9:
    sys.exit("Usage: python NN_model_prep.py parse-train.tsv cluster-train.tsv gtruth_400-train.bed exclude_region-train.bed "
             "parse-test.tsv cluster-test.tsv gtruth_400-test.bed exclude_region-test.bed")


# Obtain arguments
parse1 = argv[1]
cluster1 = argv[2]
gtruth1 = argv[3]
omit_region1 = argv[4]
parse2 = argv[5]
cluster2 = argv[6]
gtruth2 = argv[7]
omit_region2 = argv[8]


# Obtain true reads and N-gap reads for omitting
def parse_bed(cluster, gtruth, omit):
    bed = []
    with open(cluster) as f:
        for line in f:
            if len(line.split('\t')[6].split('~')) == 2:
                if len(line.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 1:
                    chrm = line.split('\t')[6].split('~')[1].split(':')[0]
                    coord1 = line.split('\t')[6].split('~')[1].split(':')[1]
                    name = line.split('\t')[8]
                    bed.append(chrm + '\t' + coord1 + '\t' + str(int(coord1) + 1) + '\t' + name)
                elif len(line.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 2:
                    chrm = line.split('\t')[6].split('~')[1].split(':')[0]
                    coord1 = line.split('\t')[6].split('~')[1].split(':')[1].split('-')[0]
                    coord2 = line.split('\t')[6].split('~')[1].split(':')[1].split('-')[1]
                    name = line.split('\t')[8]
                    bed.append(chrm + '\t' + coord1 + '\t' + str(int(coord1) + 1) + '\t' + name)
                    bed.append(chrm + '\t' + coord2 + '\t' + str(int(coord2) + 1) + '\t' + name)
            elif len(line.split('\t')[6].split('~')) == 3:
                chrm1 = line.split('\t')[6].split('~')[1].split(':')[0]
                chrm2 = line.split('\t')[6].split('~')[2].split(':')[0]
                coord1 = line.split('\t')[6].split('~')[1].split(':')[1]
                coord2 = line.split('\t')[6].split('~')[2].split(':')[1]
                name = line.split('\t')[8]
                bed.append(chrm1 + '\t' + coord1 + '\t' + str(int(coord1) + 1) + '\t' + name)
                bed.append(chrm2 + '\t' + coord2 + '\t' + str(int(coord2) + 1) + '\t' + name)
    bed_str = '\n'.join(bed)
    sv_bed = BedTool(bed_str, from_string=True)
    sv_bed = sv_bed.sort()
    truth_bed = BedTool(gtruth)
    truth_bed = truth_bed.sort()
    intersect1 = sv_bed.intersect(truth_bed, wa=True, wb=True)
    if str(omit) != 'None':
        omit_bed = BedTool(omit)
        omit_bed = omit_bed.sort()
        intersect2 = sv_bed.intersect(omit_bed, wa=True, wb=True)
    else:
        intersect2 = []
    truth_dict = {}
    omit_dict = {}
    for line in intersect1:
        line = list(line)
        truth_dict[line[3]] = 1
    for line in intersect2:
        line = list(line)
        omit_dict[line[3]] = 1
    return truth_dict, omit_dict


# For each breakpoint line, assign the long read coverage and normal cov ratio
def getoverlap(cluster, omit_dict, truth_dict):
    overlapdict = {}
    normalcovratio = {}
    classdict = {}
    with open(cluster) as f:
        for line in f:
            if line.split('\t')[8] not in omit_dict:
                overlapdict[line.split('\t')[8]] = float(line.split('\t')[10])
                normalcovratio[line.split('\t')[8]] = round(float(line.split('\t')[12])/(int(line.split('\t')[10]) +
                                                                                         float(line.split('\t')[12])), 5)
                classdict[line.split('\t')[8]] = [line.split('\t')[3].split(' ')[0], line.split('\t')[3].split(' ')[1].split('~')[0]]
                if line.split('\t')[11] != '.':
                    for d in line.split('\t')[11].split(','):
                        overlapdict[d] = float(line.split('\t')[10])
                        normalcovratio[d] = float(line.split('\t')[12])/(int(line.split('\t')[10]) + float(line.split('\t')[12]))
                        classdict[d] = [line.split('\t')[3].split(' ')[0], line.split('\t')[3].split(' ')[1].split('~')[0]]
                    if line.split('\t')[8] in truth_dict:
                        for d in line.split('\t')[11].split(','):
                            truth_dict[d] = 1
    return overlapdict, normalcovratio, classdict, truth_dict


# Filter parse data for >0 long read coverage
def filterparse(parse, overlapdict):
    newparse = []
    with open(parse) as f:
        for line in f:
            if line.split('\t')[8] in overlapdict:
                newparse.append(line)
    return newparse


# SV class detector
def svclass(sv, size):
    if str(sv) == 'S-Nov_Ins_bp' or sv == 'E-Nov_Ins_bp' or sv == 'Nov_Ins' or sv == 'Del':
        if float(size) <= 50:
            return 0
        elif float(size) <= 100:
            return 0.2
        elif float(size) <= 200:
            return 0.4
        else:
            return 1
    else:
        return 1


# Collect features for each breakpoint line and rescale the features between 0 and 1
def scalefeature(newparse, overlapdict, normalcovratio, classdict):
    signdict = {}
    for i in newparse:
        signdict[i.split('\t')[8]] = []
        for s in i.split('\t')[9].split(','):
            signdict[i.split('\t')[8]].append(s.strip('(%)'))
        signdict[i.split('\t')[8]].append(overlapdict[i.split('\t')[8]])
        signdict[i.split('\t')[8]].append(normalcovratio[i.split('\t')[8]])
        signdict[i.split('\t')[8]].append(svclass(classdict[i.split('\t')[8]][0], classdict[i.split('\t')[8]][1]))
    items = list(signdict.items())
    normsigndict = {}
    new_tup = list(map(signnorm, items))
    for i in new_tup:
        normsigndict[i[0]] = i[1]
    return normsigndict


# Set upper limits and scale signature
def signnorm(tup):
    s = tup[1]
    s[0] = round(float(s[0]) / 100, 5)
    s[1] = round(float(s[1]) / 100, 5)
    s[2] = round(float(s[2]) / 100, 5)
    s[3] = round(float(s[3]) / 100, 5)
    s[4] = round(float(s[4]) / 100, 5)
    s[5] = min(float(s[5]), 1)
    s[6] = min(float(s[6]), 1)
    s[7] = round(min(float(s[7]), 2) / 2, 5)
    s[8] = round(min(float(s[8]), 2) / 2, 5)
    s[9] = round(min(float(s[9]), 4) / 4, 5)
    s[10] = round(min(float(s[10]), 18) / 18, 5)
    s[11] = round(min(float(s[11]), 6) / 6, 5)
    s[12] = round(min(float(s[12]), 6) / 6, 5)
    s[13] = round(min(float(s[13]), 8) / 8, 5)
    s[14] = round(float(s[14]) / 100, 5)
    s[15] = round(float(s[15]) / 100, 5)
    s[16] = round(min(float(s[16]), 0.5) / 0.5, 5)
    s[17] = round(min(float(s[17]), 0.5) / 0.5, 5)
    s[18] = round(min(float(s[18]), 0.5) / 0.5, 5)
    s[19] = round(min(float(s[19]), 0.5) / 0.5, 5)
    s[20] = round(min(float(s[20]), 10) / 10, 5)
    return tup[0], s


# Label truth value to each breakpoints
def labelreads(newparse, truth_dict):
    # Add false pos reads to truth_dict
    for line in newparse:
        if line.split('\t')[8] not in truth_dict:
            truth_dict[line.split('\t')[8]] = 0
    return truth_dict


# Add truth labels to sign dict
def addlabels(normsigndict, truth_dict):
    for key in normsigndict:
        if truth_dict[key] == 1:
            normsigndict[key].append(1)
        elif truth_dict[key] == 0:
            normsigndict[key].append(0)
    return normsigndict


# Separate true and false entries and write to file
def separatenwrite(normsigndict_train, normsigndict_test):
    tmp1 = []
    tmp0 = []
    for key in normsigndict_train:
        if normsigndict_train[key][-1] == 1:
            tmp1.append(normsigndict_train[key])
        elif normsigndict_train[key][-1] == 0:
            tmp0.append(normsigndict_train[key])
    with open('train1_set.txt', 'w') as file:
        for entry in tmp1:
            file.write(','.join(str(x) for x in entry) + '\n')
    with open('train0_set.txt', 'w') as file:
        for entry in tmp0:
            file.write(','.join(str(x) for x in entry) + '\n')
    tmpts = []
    for key in normsigndict_test:
        tmpts.append(normsigndict_test[key])
    with open('test_set.txt', 'w') as file:
        for entry in tmpts:
            file.write(','.join(str(x) for x in entry) + '\n')


if __name__ == '__main__':
    truth_dict1, omit_dict1 = parse_bed(cluster1, gtruth1, omit_region1)
    truth_dict2, omit_dict2 = parse_bed(cluster2, gtruth2, omit_region2)
    overlapdict1, normalcovratio1, classdict1, truth_dict1 = getoverlap(cluster1, omit_dict1, truth_dict1)
    overlapdict2, normalcovratio2, classdict2, truth_dict2 = getoverlap(cluster2, omit_dict2, truth_dict2)
    newparse1 = filterparse(parse1, overlapdict1)
    newparse2 = filterparse(parse2, overlapdict2)
    normsigndict1 = scalefeature(newparse1, overlapdict1, normalcovratio1, classdict1)
    normsigndict2 = scalefeature(newparse2, overlapdict2, normalcovratio2, classdict2)
    truth_dict1 = labelreads(newparse1, truth_dict1)
    truth_dict2 = labelreads(newparse2, truth_dict2)
    normsigndict1 = addlabels(normsigndict1, truth_dict1)
    normsigndict2 = addlabels(normsigndict2, truth_dict2)
    separatenwrite(normsigndict1, normsigndict2)

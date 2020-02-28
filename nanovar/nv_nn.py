"""
Functions for NN model inferencing using Keras module.

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


import os
import numpy as np
from tensorflow.keras.models import load_model


# Neural network inferencing
def inference(cluster, parse, model):
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    odict, normalcovratiodict, classdict, normalcov, nmapdict, nchrdict = getoverlap(cluster)
    newpdata = filterparse(parse, odict)
    nsigndict = scalefeature(newpdata, odict, normalcovratiodict, classdict)
    tmp = []
    tmpkey = []
    for key in nsigndict:
        tmp.append(nsigndict[key])
        tmpkey.append(key)
    readarray = np.array(tmp, dtype=np.float64)
    nnmodel = load_model(model, compile=False)
    predictions = nnmodel.predict(readarray)
    predlist = [float(x[0]) for x in predictions]
    probdict = {}
    for i in range(len(tmpkey)):
        probdict[tmpkey[i]] = predlist[i]
    nn_out = []
    for i in cluster:
        finalp = finalprob(i.split('\t')[8], i.split('\t')[11], probdict)
        nn_out.append(i + '\t' + str(np.tanh(0.4 * float(odict[i.split('\t')[8]])) * finalp * readratiofilter(
            normalcovratiodict[i.split('\t')[8]], normalcov[i.split('\t')[8]]) * nmapchrfilter(nmapdict[i.split('\t')[8]])))
        # tanh(0.4B)*P where B = no. of breakend supporting reads and P = average confidence probability of breakend
        # This function alters output probability based on number of breakend supporting reads.
        # Visualize on fooplot.com
    nn_out_sort = sorted(nn_out, key=get_read_name)
    return nn_out_sort


# Additional read ratio filter post NN inference
def readratiofilter(ratio, ncov):
    if 0.65 <= ratio:
        if ncov < 6:
            score = np.exp(10 * (1 - ratio) - 3.5)
        else:
            score = np.exp(10 * (1 - ratio) - 4.5)
        return min(score, 1.0)
    else:
        return 1.0


# Additional multimapping and contig filter post NN inference
def nmapchrfilter(nmap):
    mapscore = max(-np.exp(0.3 * nmap - 5) + 1, 0)
    return mapscore


# Function for sorting list according to read name
def get_read_name(x):
    return x.split('\t')[8]


# For each breakpoint line, assign the long read coverage and normal cov ratio
def getoverlap(overlap_data):
    overlapdict = {}
    normalcovratio = {}
    classdict = {}
    normalcov = {}
    nchrdict = {}
    nmapdict = {}
    for i in overlap_data:
        overlapdict[i.split('\t')[8]] = float(i.split('\t')[10])
        normalcovratio[i.split('\t')[8]] = round(float(i.split('\t')[12])/(int(i.split('\t')[10]) + float(i.split('\t')[12])), 5)
        normalcov[i.split('\t')[8]] = float(i.split('\t')[12])
        classdict[i.split('\t')[8]] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1].split('~')[0]]
        nmapdict[i.split('\t')[8]] = int(i.split('\t')[9].split(',')[10])
        nchrdict[i.split('\t')[8]] = int(i.split('\t')[9].split(',')[12])
        if i.split('\t')[11] != '.':
            for d in i.split('\t')[11].split(','):
                overlapdict[d] = float(i.split('\t')[10])
                normalcovratio[d] = float(i.split('\t')[12])/(int(i.split('\t')[10]) + float(i.split('\t')[12]))
                normalcov[d] = float(i.split('\t')[12])
                classdict[d] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1].split('~')[0]]
    return overlapdict, normalcovratio, classdict, normalcov, nmapdict, nchrdict


# Filter parse data for >0 long read coverage
def filterparse(parse_data, overlapdict):
    newlist = []
    for i in parse_data:
        try:
            if float(overlapdict[i.split('\t')[8]]) > 0:
                newlist.append(i)
        except KeyError:
            pass
    return newlist


# SV class detector (Not in use)
# def svclass_new(sv, size):
#     if sv in ['S-Nov_Ins_bp', 'E-Nov_Ins_bp', 'Nov_Ins']:
#         if float(size) <= 50:
#             return [0.1, 0, 0, 0, 0, 0]
#         elif float(size) <= 100:
#             return [0.2, 0, 0, 0, 0, 0]
#         elif float(size) <= 200:
#             return [0.4, 0, 0, 0, 0, 0]
#         elif float(size) <= 400:
#             return [0.6, 0, 0, 0, 0, 0]
#         elif float(size) <= 1000:
#             return [0.8, 0, 0, 0, 0, 0]
#         else:
#             return [1, 0, 0, 0, 0, 0]
#     elif sv == 'Del':
#         if float(size) <= 50:
#             return [0, 0.1, 0, 0, 0, 0]
#         elif float(size) <= 100:
#             return [0, 0.2, 0, 0, 0, 0]
#         elif float(size) <= 200:
#             return [0, 0.4, 0, 0, 0, 0]
#         elif float(size) <= 400:
#             return [0, 0.6, 0, 0, 0, 0]
#         elif float(size) <= 1000:
#             return [0, 0.8, 0, 0, 0, 0]
#         else:
#             return [0, 1, 0, 0, 0, 0]
#     elif sv in ['Inv', 'Inv(1)', 'Inv(2)']:
#         return [0, 0, 1, 0, 0, 0]
#     elif sv == 'TDupl':
#         return [0, 0, 0, 1, 0, 0]
#     elif sv in ['Intra-Ins', 'Intra-Ins(1)', 'Intra-Ins(2)']:
#         return [0, 0, 0, 0, 1, 0]
#     elif sv in ['InterTx', 'Inter-Ins(1)', 'Inter-Ins(2)']:
#         return [0, 0, 0, 0, 0, 1]
#     else:
#         raise Exception('ERROR: SV class not recognised %s' % sv)


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
def scalefeature(parse_data, overlapdict, normalcovratio, classdict):
    signdict = {}
    for i in parse_data:
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
    s[10] = round((min(float(s[10]), 18) - 1) / 18, 5)
    s[11] = round((min(float(s[11]), 6) - 1) / 6, 5)
    s[12] = round((min(float(s[12]), 4) - 1) / 4, 5)
    s[13] = round(min(float(s[13]), 8) / 8, 5)
    s[14] = round(float(s[14]) / 100, 5)
    s[15] = round(float(s[15]) / 100, 5)
    s[16] = round(min(float(s[16]), 0.5) / 0.5, 5)
    s[17] = round(min(float(s[17]), 0.5) / 0.5, 5)
    s[18] = round(min(float(s[18]), 0.5) / 0.5, 5)
    s[19] = round(min(float(s[19]), 0.5) / 0.5, 5)
    s[20] = round(min(float(s[20]), 10) / 10, 5)
    return tup[0], s


# Function to average out NN score
def finalprob(leadname, fellow, probdict):
    leadprob = probdict[leadname]
    prob = float(leadprob)
    if fellow != '.':
        fellows = fellow.split(',')
        n = len(fellows) + 1
        for i in fellows:
            prob += float(probdict[i])
    elif fellow == '.':
        n = 1
    else:
        raise Exception("Error: Breakend fellow symbol unrecognised")
    fprob = float(prob)/n
    return fprob


# Output SV read overlap partners
def svread_ovl(wk_dir, out_nn):
    data = open(os.path.join(wk_dir, 'sv_support_reads.tsv'), 'w')
    data.write('SV-ID\tSupporting_reads (readname~index1,readname~index2...)\n')
    for line in out_nn:
        svid = line.split('\t')[6].split('~')[0]
        data.write(svid + '\t' + line.split('\t')[8] + ',' + line.split('\t')[11] + '\n')
    data.close()

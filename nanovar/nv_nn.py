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
import sys
import numpy as np
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
from keras.models import load_model
sys.stderr = stderr


# Neural network inferencing
def inference(cluster, parse, model):
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    odict, normalcovratiodict, classdict = getoverlap(cluster)
    newpdata = filterparse(parse, odict)
    nsigndict = scalefeature(newpdata, odict, normalcovratiodict, classdict)
    tmp = []
    tmpkey = []
    for key in nsigndict:
        tmp.append(nsigndict[key])
        tmpkey.append(key)
    readarray = np.array(tmp, dtype=np.float64)
    nnmodel = load_model(model)
    predictions = nnmodel.predict(readarray)
    predlist = [float(x[0]) for x in predictions]
    probdict = {}
    for i in range(len(tmpkey)):
        probdict[tmpkey[i]] = predlist[i]
    nn_out = []
    for i in cluster:
        finalp = finalprob(i.split('\t')[8], i.split('\t')[11], probdict)
        nn_out.append(i + '\t' + str(np.tanh(0.4 * float(odict[i.split('\t')[8]])) * finalp))
        # tanh(0.4B)*P where B = no. of breakend supporting reads and P = average confidence probability of breakend
        # This function alters output probability based on number of breakend supporting reads.
        # Visualize on fooplot.com
    nn_out_sort = sorted(nn_out, key=get_read_name)
    return nn_out_sort


# Function for sorting list according to read name
def get_read_name(x):
    return x.split('\t')[8]


# For each breakpoint line, assign the long read coverage and normal cov ratio
def getoverlap(overlap_data):
    overlapdict = {}
    normalcovratio = {}
    classdict = {}
    for i in overlap_data:
        overlapdict[i.split('\t')[8]] = i.split('\t')[10]
        normalcovratio[i.split('\t')[8]] = float(i.split('\t')[12])/(int(i.split('\t')[10]) + float(i.split('\t')[12]))
        classdict[i.split('\t')[8]] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1].split('~')[0]]
        if i.split('\t')[11] != '.':
            for d in i.split('\t')[11].split(','):
                overlapdict[d] = i.split('\t')[10]
                normalcovratio[d] = float(i.split('\t')[12])/(int(i.split('\t')[10]) + float(i.split('\t')[12]))
                classdict[d] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1].split('~')[0]]
    return overlapdict, normalcovratio, classdict


# Filter parse data for >0 long read coverage
def filterparse(parse_data, overlapdict):
    newlist = []
    for i in parse_data:
        try:
            if int(overlapdict[i.split('\t')[8]]) > 0:
                newlist.append(i)
        except KeyError:
            pass
    return newlist


# Del and Ins detector
def delins(sv):
    if str(sv) == 'S-Nov_Ins_bp' or sv == 'E-Nov_Ins_bp' or sv == 'Nov_Ins' or sv == 'Del':
        return True


# Collect features for each breakpoint line and rescale the features between 0 and 1
def scalefeature(parse_data, overlapdict, normalcovratio, classdict):
    signdict = {}
    for i in parse_data:
        signdict[i.split('\t')[8]] = []
        for s in i.split('\t')[9].split(','):
            signdict[i.split('\t')[8]].append(s.strip('(%)'))
        signdict[i.split('\t')[8]].append(overlapdict[i.split('\t')[8]])
        signdict[i.split('\t')[8]].append(normalcovratio[i.split('\t')[8]])
        if delins(classdict[i.split('\t')[8]][0]):
            if float(classdict[i.split('\t')[8]][1]) <= 50:
                signdict[i.split('\t')[8]].append(0)
            elif float(classdict[i.split('\t')[8]][1]) <= 100:
                signdict[i.split('\t')[8]].append(0.2)
            elif float(classdict[i.split('\t')[8]][1]) <= 200:
                signdict[i.split('\t')[8]].append(0.4)
            else:
                signdict[i.split('\t')[8]].append(1)
        else:
            signdict[i.split('\t')[8]].append(1)
    le = len(signdict[parse_data[0].split('\t')[8]])
    minmaxlist = []
    for i in range(le):
        tmp = []
        for key in signdict:
            tmp.append(float(signdict[key][i]))
        minmaxlist.append([min(tmp), max(tmp)])
    normsigndict = {}
    for key in signdict:
        normsigndict[key] = []
        for i in range(le):
            normsigndict[key].append(str((float(signdict[key][i]) - minmaxlist[i][0])/(minmaxlist[i][1] - minmaxlist[i][0]
                                                                                       + 0.0001)))
    return normsigndict


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
    data = open(os.path.join(wk_dir, 'svread-overlap.tsv'), 'w')
    data.write('Leader_index~readname\tSupporting_readname~index1,Supporting_readname~index2...')
    for line in out_nn:
        lead = line.split('\t')[8].split('~')[1] + '~' + line.split('\t')[8].split('~')[0]
        data.write(lead + '\t' + line.split('\t')[11] + '\n')
    data.close()

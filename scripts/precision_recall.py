
### Calculate precision and recall and plot PRC curve of homozygous SVs for each tool

from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib
import numpy as np

# Generate PRC with scoring
def prc_with_score(callset, tp, fn):
    fn_data = open(fn, 'r').read().splitlines()
    # Make readname to bp dict and breakpoint_name to score dict
    rnamebpdict = {}
    bpscoredict = {}
    with open(tp) as tp_data:
        for line in tp_data:
            rname = line.split('\t')[0].strip()
            score = line.split('\t')[1].strip()
            bpname = line.split('\t')[2].strip()
            try:
                rnamebpdict[rname].append(bpname)
            except:
                rnamebpdict[rname] = []
                rnamebpdict[rname].append(bpname)
            try:
                bpscoredict[bpname].append(float(score))
            except:
                bpscoredict[bpname] = []
                bpscoredict[bpname].append(float(score))
    #Make scorelist for scaling and novodict 
    scorelist = []
    sdict = {}
    for key in bpscoredict:
        scorelist.append(float(max(bpscoredict[key])))
    with open(callset) as cs:
        for line in cs:
            rname = line.split('\t')[0].strip()
            score = line.split('\t')[1].strip()
            try:
                if rnamebpdict[rname]:
                    pass
            except KeyError:
                scorelist.append(float(score))
    minmaxscore = [min(scorelist), max(scorelist)]
    # Make dict bpname/readname to score and truth value
    sdict = {}
    with open(callset) as cs:
        for line in cs:
            rname = line.split('\t')[0].strip()
            score = line.split('\t')[1].strip()
            try:
                if rnamebpdict[rname]:
                    for i in rnamebpdict[rname]:
                        sdict[i] = []
                        sdict[i].append((max(bpscoredict[i]) - minmaxscore[0])/(minmaxscore[1]+0.01 - minmaxscore[0]))
                        sdict[i].append(float(1))
            except KeyError:
                sdict[rname] = []
                sdict[rname].append((float(score) - minmaxscore[0])/(minmaxscore[1]+0.01 - minmaxscore[0]))
                sdict[rname].append(float(0))
    tmp = []
    for key in sdict:
        tmp.append(sdict[key])
    readarray = np.array(tmp, dtype=np.float64)
    pred = np.split(readarray, [1], axis=1)[0]
    actual = np.split(readarray, [1], axis=1)[1]
    for line in fn_data:
        pred = np.append(pred, [[-1]], axis=0)
        actual = np.append(actual, [[1]], axis=0)
    precision, recall, threshold = precision_recall_curve(actual, pred) # Generate Precision recall values
    np.put(precision, [0], precision[1]) # Replaces first value with second value for asethetics
    np.put(recall, [0], recall[1]) # Replaces first value with second value for asethetics
    auc_value = auc(recall, precision) # Get area under curve value
    return precision, recall, threshold, auc_value

# Generate PRC without scoring
def prc_no_score(callset, tp, fn):
    # Read false neg data
    fn_data = open(fn, 'r').read().splitlines()
    # Make readname to bp dict
    rnamebpdict = {}
    with open(tp) as tp_data:
        for line in tp_data:
            rname = line.split('\t')[0].strip()
            bpname = line.split('\t')[2].strip()
            try:
                rnamebpdict[rname].append(bpname)
            except:
                rnamebpdict[rname] = []
                rnamebpdict[rname].append(bpname)
    # Make dict bpname/readname to score and truth value
    sdict = {}
    with open(callset) as cs:
        for line in cs:
            rname = line.split('\t')[0].strip()
            try:
                if rnamebpdict[rname]:
                    for i in rnamebpdict[rname]:
                        sdict[i] = []
                        sdict[i].append(float(1))
                        sdict[i].append(float(1))
            except KeyError:
                sdict[rname] = []
                sdict[rname].append(float(1))
                sdict[rname].append(float(0))
    tmp = []
    for key in sdict:
        tmp.append(sdict[key])
    readarray = np.array(tmp, dtype=np.float64)
    pred = np.split(readarray, [1], axis=1)[0]
    actual = np.split(readarray, [1], axis=1)[1]
    for line in fn_data:
        pred = np.append(pred, [[-1]], axis=0)
        actual = np.append(actual, [[1]], axis=0)
    precision, recall, threshold = precision_recall_curve(actual, pred) # Generate Precision recall values
    np.put(precision, [0], [0]) # Replaces first value with 0 for asethetics
    np.put(recall, [0], recall[1]) # Replaces first value with second value for asethetics
    precision = precision[1]
    recall = recall[1]
    return precision, recall, threshold

# NanoVar (with thresholds)
callset = './combine1/nanovar/total.callset'
tp = './combine1/nanovar/total-tp'
fn = './combine1/nanovar/total-fn'

precision_nv1, recall_nv1, threshold_nv1, auc_nv1 = prc_with_score(callset, tp, fn)

callset = './combine2/nanovar/total.callset'
tp = './combine2/nanovar/total-tp'
fn = './combine2/nanovar/total-fn'

precision_nv2, recall_nv2, threshold_nv2, auc_nv2 = prc_with_score(callset, tp, fn)

callset = './combine3/nanovar/total.callset'
tp = './combine3/nanovar/total-tp'
fn = './combine3/nanovar/total-fn'

precision_nv3, recall_nv3, threshold_nv3, auc_nv3 = prc_with_score(callset, tp, fn)

ave_auc_nv = (auc_nv1 + auc_nv2 + auc_nv3)/3

# threshold = 1.0
precision_nv_t1 = precision_nv1[11]
recall_nv_t1 = recall_nv1[11]

precision_nv_t2 = precision_nv2[11]
recall_nv_t2 = recall_nv2[11]

precision_nv_t3 = precision_nv3[11]
recall_nv_t3 = recall_nv3[11]



# Picky (without thresholds)
callset = './combine1/picky/total.callset'
tp = './combine1/picky/total-tp'
fn = './combine1/picky/total-fn'

precision_p1, recall_p1, threshold_p1 = prc_no_score(callset, tp, fn)

callset = './combine2/picky/total.callset'
tp = './combine2/picky/total-tp'
fn = './combine2/picky/total-fn'

precision_p2, recall_p2, threshold_p2 = prc_no_score(callset, tp, fn)

callset = './combine3/picky/total.callset'
tp = './combine3/picky/total-tp'
fn = './combine3/picky/total-fn'

precision_p3, recall_p3, threshold_p3 = prc_no_score(callset, tp, fn)



# Sniffles (without thresholds)
callset = './combine1/sniffles/total.callset'
tp = './combine1/sniffles/total-tp'
fn = './combine1/sniffles/total-fn'

precision_snif1, recall_snif1, threshold_snif1 = prc_no_score(callset, tp, fn)

callset = './combine2/sniffles/total.callset'
tp = './combine2/sniffles/total-tp'
fn = './combine2/sniffles/total-fn'

precision_snif2, recall_snif2, threshold_snif2 = prc_no_score(callset, tp, fn)

callset = './combine3/sniffles/total.callset'
tp = './combine3/sniffles/total-tp'
fn = './combine3/sniffles/total-fn'

precision_snif3, recall_snif3, threshold_snif3 = prc_no_score(callset, tp, fn)



# NanoSV (with thresholds)
callset = './combine1/nanoSV/total.callset'
tp = './combine1/nanoSV/total-tp'
fn = './combine1/nanoSV/total-fn'

precision_ns1, recall_ns1, threshold_ns1, auc_ns1 = prc_with_score(callset, tp, fn)

callset = './combine2/nanoSV/total.callset'
tp = './combine2/nanoSV/total-tp'
fn = './combine2/nanoSV/total-fn'

precision_ns2, recall_ns2, threshold_ns2, auc_ns2 = prc_with_score(callset, tp, fn)

callset = './combine3/nanoSV/total.callset'
tp = './combine3/nanoSV/total-tp'
fn = './combine3/nanoSV/total-fn'

precision_ns3, recall_ns3, threshold_ns3, auc_ns3 = prc_with_score(callset, tp, fn)

ave_auc_ns = (auc_ns1 + auc_ns2 + auc_ns3)/3

# threshold = 0
precision_ns_t1 = precision_ns1[0]
recall_ns_t1 = recall_ns1[0]

precision_ns_t2 = precision_ns2[0]
recall_ns_t2 = recall_ns2[0]

precision_ns_t3 = precision_ns3[0]
recall_ns_t3 = recall_ns3[0]



# SVIM (with thresholds)
callset = './combine1/svim/total.callset'
tp = './combine1/svim/total-tp'
fn = './combine1/svim/total-fn'

precision_sm1, recall_sm1, threshold_sm1, auc_sm1 = prc_with_score(callset, tp, fn)

callset = './combine2/svim/total.callset'
tp = './combine2/svim/total-tp'
fn = './combine2/svim/total-fn'

precision_sm2, recall_sm2, threshold_sm2, auc_sm2 = prc_with_score(callset, tp, fn)

callset = './combine3/svim/total.callset'
tp = './combine3/svim/total-tp'
fn = './combine3/svim/total-fn'

precision_sm3, recall_sm3, threshold_sm3, auc_sm3 = prc_with_score(callset, tp, fn)

ave_auc_sm = (auc_sm1 + auc_sm2 + auc_sm3)/3

# threshold = 0
precision_sm_t1 = precision_sm1[0]
recall_sm_t1 = recall_sm1[0]

precision_sm_t2 = precision_sm2[0]
recall_sm_t2 = recall_sm2[0]

precision_sm_t3 = precision_sm3[0]
recall_sm_t3 = recall_sm3[0]



# novobreak (with thresholds)
callset = './combine1/novobreak/total.callset'
tp = './combine1/novobreak/total-tp'
fn = './combine1/novobreak/total-fn'

precision_nb1, recall_nb1, threshold_nb1, auc_nb1 = prc_with_score(callset, tp, fn)

callset = './combine2/novobreak/total.callset'
tp = './combine2/novobreak/total-tp'
fn = './combine2/novobreak/total-fn'

precision_nb2, recall_nb2, threshold_nb2, auc_nb2 = prc_with_score(callset, tp, fn)

callset = './combine3/novobreak/total.callset'
tp = './combine3/novobreak/total-tp'
fn = './combine3/novobreak/total-fn'

precision_nb3, recall_nb3, threshold_nb3, auc_nb3 = prc_with_score(callset, tp, fn)

ave_auc_nb = (auc_nb1 + auc_nb2 + auc_nb3)/3

# threshold = 27.5
precision_nb_t1 = precision_nb1[56]
recall_nb_t1 = recall_nb1[56]

precision_nb_t2 = precision_nb2[56]
recall_nb_t2 = recall_nb2[56]

precision_nb_t3 = precision_nb3[56]
recall_nb_t3 = recall_nb3[56]



# Delly (without threshold)
callset = './combine1/delly/total.callset'
tp = './combine1/delly/total-tp'
fn = './combine1/delly/total-fn'

precision_de1, recall_de1, threshold_de1 = prc_no_score(callset, tp, fn)

callset = './combine2/delly/total.callset'
tp = './combine2/delly/total-tp'
fn = './combine2/delly/total-fn'

precision_de2, recall_de2, threshold_de2 = prc_no_score(callset, tp, fn)

callset = './combine3/delly/total.callset'
tp = './combine3/delly/total-tp'
fn = './combine3/delly/total-fn'

precision_de3, recall_de3, threshold_de3 = prc_no_score(callset, tp, fn)



matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
params = {'axes.labelsize': 9,'axes.titlesize':8, 'legend.fontsize': 7, 'xtick.labelsize':8, 'ytick.labelsize': 8}
matplotlib.rcParams.update(params)



# Plot PRC
fig = plt.figure(figsize=(3.7, 3.9))
ax = fig.add_subplot(111)

plt.plot(recall_nv_t1, precision_nv_t1, color='mediumslateblue', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_nv_t2, precision_nv_t2, color='mediumslateblue', alpha=1.0, marker='v', markersize=4, linestyle='None')
nv, = plt.plot(recall_nv_t3, precision_nv_t3, label="NanoVar", color='mediumslateblue', alpha=1.0, marker='v', markersize=4, linestyle='None')

plt.plot(recall_sm_t1, precision_sm_t1, color='saddlebrown', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_sm_t2, precision_sm_t2, color='saddlebrown', alpha=1.0, marker='v', markersize=4, linestyle='None')
sm, = plt.plot(recall_sm_t3, precision_sm_t3, label="SVIM", color='saddlebrown', alpha=1.0, marker='v', markersize=4, linestyle='None')

plt.plot(recall_ns_t1, precision_ns_t1, color='#ffd919', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_ns_t2, precision_ns_t2, color='#ffd919', alpha=1.0, marker='v', markersize=4, linestyle='None')
ns, = plt.plot(recall_ns_t3, precision_ns_t3, label="NanoSV", color='#ffd919', alpha=1.0, marker='v', markersize=4, linestyle='None')

plt.plot(recall_nb_t1, precision_nb_t1, color='crimson', alpha=1.0, marker='o', markersize=4, linestyle='None')
plt.plot(recall_nb_t2, precision_nb_t2, color='crimson', alpha=1.0, marker='o', markersize=4, linestyle='None')
nb, = plt.plot(recall_nb_t3, precision_nb_t3, label="novoBreak", color='crimson', alpha=1.0, marker='o', markersize=4, linestyle='None')

pick, = plt.plot(recall_p1, precision_p1, label="Picky", color='#3ca63b', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_p2, precision_p2, color='#3ca63b', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_p3, precision_p3, color='#3ca63b', alpha=1.0, marker='v', markersize=4, linestyle='None')

snif, = plt.plot(recall_snif1, precision_snif1, label="Sniffles", color='#23a7e4', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_snif2, precision_snif2, color='#23a7e4', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_snif3, precision_snif3, color='#23a7e4', alpha=1.0, marker='v', markersize=4, linestyle='None')

dell, = plt.plot(recall_de1, precision_de1, label="Delly", color='darkorange', alpha=1.0, marker='o', markersize=4, linestyle='None')
plt.plot(recall_de2, precision_de2, color='darkorange', alpha=1.0, marker='o', markersize=4, linestyle='None')
plt.plot(recall_de3, precision_de3, color='darkorange', alpha=1.0, marker='o', markersize=4, linestyle='None')

# Markers for legend
marker_4x = mlines.Line2D([], [], color='black', marker='v', linestyle='None', markersize=4, label='4X depth')
marker_8x = mlines.Line2D([], [], color='black', marker='D', linestyle='None', markersize=4, label='8X depth')
marker_12x = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=4, label='12X depth')

# F1 score curves
f_scores = np.linspace(0.2, 0.8, num=4)
f_scores = np.append(f_scores,0.9)
lines = []
labels = []
for f_score in f_scores:
    x = np.linspace(0.001, 1)
    y = f_score * x / (2 * x - f_score)
    fline, = plt.plot(x[0<=y], y[0<=y], label="$F_1$ score curves", color='gray', alpha=0.4)
    plt.annotate('$f_1$={0:0.1f}'.format(f_score), xy=(0.93 - 0.01, y[46] + 0.02), fontsize=5.5)

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.yticks(np.arange(0.0, 1.05, 0.1))
plt.xticks(np.arange(0.0, 1.05, 0.1))
fline.set_alpha(0.8)
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.95])

# Create legends
legend1 = plt.legend(handles=[nv, pick, ns, sm, snif, marker_4x, marker_8x, marker_12x], bbox_to_anchor=(0.45,-0.215), ncol=2, fancybox=True)
legend1.get_frame().set_edgecolor('black')
ax = plt.gca().add_artist(legend1)
legend2 = plt.legend(handles=[nb, dell], bbox_to_anchor=(1.05,-0.215), ncol=2, fancybox=True)
legend2.get_frame().set_edgecolor('black')
ax = plt.gca().add_artist(legend2)
plt.legend(handles=[fline], bbox_to_anchor=(0.97,-0.293), frameon=False)
plt.text(-0.003,-0.228, "Long-read SV callers", fontsize=8)
plt.text(0.590,-0.228, "Short-read SV callers", fontsize=8)
fig.tight_layout(rect=[0, 0.1, 1, 1.05])
plt.savefig('fig.svg', format='svg', dpi = 350)




### Calculate precision and recall and plot PRC curve of heterozygous SVs for each tool

# NanoVar
callset = './combine4x/nanovar/total.callset'
tp = './combine4x/nanovar/total-tp'
fn = './combine4x/nanovar/total-fn'

precision_nv1, recall_nv1, threshold_nv1, auc_nv1 = prc_with_score(callset, tp, fn)

callset = './combine8x/nanovar/total.callset'
tp = './combine8x/nanovar/total-tp'
fn = './combine8x/nanovar/total-fn'

precision_nv2, recall_nv2, threshold_nv2, auc_nv2 = prc_with_score(callset, tp, fn)

callset = './combine12x/nanovar/total.callset'
tp = './combine12x/nanovar/total-tp'
fn = './combine12x/nanovar/total-fn'

precision_nv3, recall_nv3, threshold_nv3, auc_nv3 = prc_with_score(callset, tp, fn)

ave_auc_nv = (auc_nv1 + auc_nv2 + auc_nv3)/3

# threshold = 1.0
precision_nv_t1 = precision_nv1[11]
recall_nv_t1 = recall_nv1[11]

precision_nv_t2 = precision_nv2[11]
recall_nv_t2 = recall_nv2[11]

precision_nv_t3 = precision_nv3[11]
recall_nv_t3 = recall_nv3[11]



# Picky (without thresholds)
callset = './combine4x/picky/total.callset'
tp = './combine4x/picky/total-tp'
fn = './combine4x/picky/total-fn'

precision_p1, recall_p1, threshold_p1 = prc_no_score(callset, tp, fn)

callset = './combine8x/picky/total.callset'
tp = './combine8x/picky/total-tp'
fn = './combine8x/picky/total-fn'

precision_p2, recall_p2, threshold_p2 = prc_no_score(callset, tp, fn)

callset = './combine12x/picky/total.callset'
tp = './combine12x/picky/total-tp'
fn = './combine12x/picky/total-fn'

precision_p3, recall_p3, threshold_p3 = prc_no_score(callset, tp, fn)



# Sniffles (without thresholds)
callset = './combine4x/sniffles/total.callset'
tp = './combine4x/sniffles/total-tp'
fn = './combine4x/sniffles/total-fn'

precision_snif1, recall_snif1, threshold_snif1 = prc_no_score(callset, tp, fn)

callset = './combine8x/sniffles/total.callset'
tp = './combine8x/sniffles/total-tp'
fn = './combine8x/sniffles/total-fn'

precision_snif2, recall_snif2, threshold_snif2 = prc_no_score(callset, tp, fn)

callset = './combine12x/sniffles/total.callset'
tp = './combine12x/sniffles/total-tp'
fn = './combine12x/sniffles/total-fn'

precision_snif3, recall_snif3, threshold_snif3 = prc_no_score(callset, tp, fn)



# NanoSV (with thresholds)
callset = './combine4x/nanoSV/total.callset'
tp = './combine4x/nanoSV/total-tp'
fn = './combine4x/nanoSV/total-fn'

precision_ns1, recall_ns1, threshold_ns1, auc_ns1 = prc_with_score(callset, tp, fn)

callset = './combine8x/nanoSV/total.callset'
tp = './combine8x/nanoSV/total-tp'
fn = './combine8x/nanoSV/total-fn'

precision_ns2, recall_ns2, threshold_ns2, auc_ns2 = prc_with_score(callset, tp, fn)

callset = './combine12x/nanoSV/total.callset'
tp = './combine12x/nanoSV/total-tp'
fn = './combine12x/nanoSV/total-fn'

precision_ns3, recall_ns3, threshold_ns3, auc_ns3 = prc_with_score(callset, tp, fn)

ave_auc_ns = (auc_ns1 + auc_ns2 + auc_ns3)/3

# threshold = 0
precision_ns_t1 = precision_ns1[0]
recall_ns_t1 = recall_ns1[0]

precision_ns_t2 = precision_ns2[0]
recall_ns_t2 = recall_ns2[0]

precision_ns_t3 = precision_ns3[0]
recall_ns_t3 = recall_ns3[0]



# SVIM (with thresholds)
callset = './combine4x/svim/total.callset'
tp = './combine4x/svim/total-tp'
fn = './combine4x/svim/total-fn'

precision_sm1, recall_sm1, threshold_sm1, auc_sm1 = prc_with_score(callset, tp, fn)

callset = './combine8x/svim/total.callset'
tp = './combine8x/svim/total-tp'
fn = './combine8x/svim/total-fn'

precision_sm2, recall_sm2, threshold_sm2, auc_sm2 = prc_with_score(callset, tp, fn)

callset = './combine12x/svim/total.callset'
tp = './combine12x/svim/total-tp'
fn = './combine12x/svim/total-fn'

precision_sm3, recall_sm3, threshold_sm3, auc_sm3 = prc_with_score(callset, tp, fn)

ave_auc_sm = (auc_sm1 + auc_sm2 + auc_sm3)/3

# threshold = 0
precision_sm_t1 = precision_sm1[0]
recall_sm_t1 = recall_sm1[0]

precision_sm_t2 = precision_sm2[0]
recall_sm_t2 = recall_sm2[0]

precision_sm_t3 = precision_sm3[0]
recall_sm_t3 = recall_sm3[0]



# novobreak (with thresholds)
callset = './combine_novobreak/total.callset'
tp = './combine_novobreak/total-tp'
fn = './combine_novobreak/total-fn'

precision_nb1, recall_nb1, threshold_nb1, auc_nb1 = prc_with_score(callset, tp, fn)

# threshold = 27.5
precision_nb_t1 = precision_nb1[56]
recall_nb_t1 = recall_nb1[56]



# Delly (without threshold)
callset = './combine_delly/total.callset'
tp = './combine_delly/total-tp'
fn = './combine_delly/total-fn'

precision_de1, recall_de1, threshold_de1 = prc_no_score(callset, tp, fn)



matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
params = {'axes.labelsize': 9,'axes.titlesize':8, 'legend.fontsize': 7, 'xtick.labelsize':8, 'ytick.labelsize': 8}
matplotlib.rcParams.update(params)


# Plot PRC
fig = plt.figure(figsize=(3.7, 3.9))
ax = fig.add_subplot(111)

nv, = plt.plot(recall_nv_t1, precision_nv_t1, label="NanoVar", color='mediumslateblue', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_nv_t2, precision_nv_t2, color='mediumslateblue', alpha=1.0, marker='D', markersize=4, linestyle='None')
plt.plot(recall_nv_t3, precision_nv_t3, color='mediumslateblue', alpha=1.0, marker='*', markersize=4, linestyle='None')

sm, = plt.plot(recall_sm_t1, precision_sm_t1, label="SVIM", color='saddlebrown', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_sm_t2, precision_sm_t2, color='saddlebrown', alpha=1.0, marker='D', markersize=4, linestyle='None')
plt.plot(recall_sm_t3, precision_sm_t3, color='saddlebrown', alpha=1.0, marker='*', markersize=4, linestyle='None')

ns, = plt.plot(recall_ns_t1, precision_ns_t1, label="NanoSV", color='#ffd919', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_ns_t2, precision_ns_t2, color='#ffd919', alpha=1.0, marker='D', markersize=4, linestyle='None')
plt.plot(recall_ns_t3, precision_ns_t3, color='#ffd919', alpha=1.0, marker='*', markersize=4, linestyle='None')

nb, = plt.plot(recall_nb_t1, precision_nb_t1, label="novoBreak (53X)", color='crimson', alpha=1.0, marker='o', markersize=4, linestyle='None')

dell, = plt.plot(recall_de1, precision_de1, label="Delly (53X)", color='darkorange', alpha=1.0, marker='o', markersize=4, linestyle='None')

pick, = plt.plot(recall_p1, precision_p1, label="Picky", color='#3ca63b', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_p2, precision_p2, color='#3ca63b', alpha=1.0, marker='D', markersize=4, linestyle='None')
plt.plot(recall_p3, precision_p3, color='#3ca63b', alpha=1.0, marker='*', markersize=4, linestyle='None')

snif, = plt.plot(recall_snif1, precision_snif1, label="Sniffles", color='#23a7e4', alpha=1.0, marker='v', markersize=4, linestyle='None')
plt.plot(recall_snif2, precision_snif2, color='#23a7e4', alpha=1.0, marker='D', markersize=4, linestyle='None')
plt.plot(recall_snif3, precision_snif3, color='#23a7e4', alpha=1.0, marker='*', markersize=4, linestyle='None')

# Markers for legend
marker_4x = mlines.Line2D([], [], color='black', marker='v', linestyle='None', markersize=4, label='4X depth')
marker_8x = mlines.Line2D([], [], color='black', marker='D', linestyle='None', markersize=4, label='8X depth')
marker_12x = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=4, label='12X depth')

# F1 score curves
f_scores = np.linspace(0.2, 0.8, num=4)
f_scores = np.append(f_scores,0.9)
lines = []
labels = []
for f_score in f_scores:
    x = np.linspace(0.001, 1)
    y = f_score * x / (2 * x - f_score)
    fline, = plt.plot(x[0<=y], y[0<=y], label="$F_1$ score curves", color='gray', alpha=0.4)
    plt.annotate('$f_1$={0:0.1f}'.format(f_score), xy=(0.93 - 0.01, y[46] + 0.02), fontsize=5.5)

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.yticks(np.arange(0.0, 1.05, 0.1))
plt.xticks(np.arange(0.0, 1.05, 0.1))
fline.set_alpha(0.8)
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.95])

legend1 = plt.legend(handles=[nv, pick, ns, sm, snif, marker_4x, marker_8x, marker_12x], bbox_to_anchor=(0.45,-0.215), ncol=2, fancybox=True)
legend1.get_frame().set_edgecolor('black')
ax = plt.gca().add_artist(legend1)
legend2 = plt.legend(handles=[nb, dell], bbox_to_anchor=(0.97,-0.215), ncol=1, fancybox=True)
legend2.get_frame().set_edgecolor('black')
ax = plt.gca().add_artist(legend2)
plt.legend(handles=[fline], bbox_to_anchor=(0.97,-0.365), frameon=False)
plt.text(-0.003,-0.228, "Long-read SV callers", fontsize=8)
plt.text(0.590,-0.228, "Short-read SV callers", fontsize=8)
fig.tight_layout(rect=[0, 0.1, 1, 1.05])
plt.savefig('fig.svg', format='svg', dpi = 350)

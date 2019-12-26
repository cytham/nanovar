'''
# Shell commands
# Parse ground-truth deletions
tail -n +2 deletions.csv | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' | bedtools sort -i - > del.bed

# Parse ground-truth inversions
tail -n +2 inversions.csv | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' | bedtools sort -i - > inv.bed

# Parse ground-truth tandem duplications
tail -n +2 tandemDuplications.csv | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' | bedtools sort -i - > dup.bed

# Parse ground-truth transpositions
tail -n +2 insertions.csv | awk '{if ($9=="FALSE") print $2"\t"$3"\t"$4"\t"$8"\t""DEL_"$1"\n"$5"\t"$6"\t"$6+1"\t"$8"\t"$1;else print $5"\t"$6"\t"$6+1"\t"$8"\t"$1}' | bedtools sort -i - > tpo.bed

# Parse ground-truth viral insertions
tail -n +2 insertions.csv | cut -f 1,5-8 | awk -F'\t' '{print $2"\t"$3"\t"$3+1"\t"$5"\tins"$1}' | bedtools sort -i - > ins.bed

# Get repeatmasker hg19 from UCSC tables and parse into a BED file
tail -n +2 hg19_repeatmasker.tsv | cut -f 6,7,8,10,11,12,13 > hg19_repeatmasker.bed

# Combine 32k SV and 10k SV into a single ground truth BED file and intersect with repeat masker BED file
cd combine1
cat ../del.bed ../inv.bed ../dup.bed ../tpo.bed ../ins.bed | bedtools sort -i - > total_sv_actual.bed
bedtools intersect -wa -wb -a total_sv_actual.bed -b ./hg19_repeatmasker.bed | cut -f 5,10-12 > total_sv_actual.names-repeats.tsv
'''


# Python commands
# Define ground truth files
gtruth1 = "total_sv_actual.bed"
gtruth_rep1 = "total_sv_actual.names-repeats.tsv"

# Define true positive files of tools for Sim1 32k + Sim1 10k homozygous SV (4X for long-read, 53X for short-read)
nv1homo = "./homo/combine1/nanovar/total-tp"
snif1homo = "./homo/combine1/sniffles/total-tp"
sv1homo = "./homo/combine1/svim/total-tp"
ns1homo = "./homo/combine1/nanoSV/total-tp"
pk1homo = "./homo/combine1/picky/total-tp"
nb1homo = "./homo/combine1/novobreak/total-tp"
de1homo = "./homo/combine1/delly/total-tp"

# Define true positive files of tools for Sim1 32k + Sim1 10k heterozygous SV (4X for long-read, 53X for short-read)
nv1het = "./hetero/combine4x/nanovar/total-tp"
snif1het = "./hetero/combine4x/sniffles/total-tp"
sv1het = "./hetero/combine4x/svim/total-tp"
ns1het = "./hetero/combine4x/nanoSV/total-tp"
pk1het = "./hetero/combine4x/picky/total-tp"
nb1het = "./hetero/combine_novobreak/total-tp"
de1het = "./hetero/combine_delly/total-tp"

# Define types of repeat families
rep_types = ["DNA", "LINE", "SINE", "LTR", "Low_complexity", "Satellite", "Simple_repeat"]

# Extract repeats information in ground-truth SVs
sv_rep_file = open(gtruth_rep1, 'r').read().splitlines()
sv_rep_dict = {}
for line in sv_rep_file:
    if line.split('\t')[2] in rep_types:
        try:
            sv_rep_dict[line.split('\t')[0]].append(line.split('\t')[2])
            sv_rep_dict[line.split('\t')[0]] = list(set(sv_rep_dict[line.split('\t')[0]]))
        except KeyError:
            sv_rep_dict[line.split('\t')[0]] = []
            sv_rep_dict[line.split('\t')[0]].append(line.split('\t')[2])
    else:
        try:
            sv_rep_dict[line.split('\t')[0]].append('others')
            sv_rep_dict[line.split('\t')[0]] = list(set(sv_rep_dict[line.split('\t')[0]]))
        except KeyError:
            sv_rep_dict[line.split('\t')[0]] = []
            sv_rep_dict[line.split('\t')[0]].append('others')

total_file = open(gtruth1, 'r').read().splitlines()
for line in total_file:
    try:
        if sv_rep_dict[line.split('\t')[4]]:
            pass
    except KeyError:
        sv_rep_dict[line.split('\t')[4]] = ['non-repeat']

# Count total number of each repeat families
total_rep = {"DNA": 0, "LINE": 0, "SINE": 0, "LTR": 0, "Low_complexity": 0, "Satellite": 0, "Simple_repeat": 0, 'others': 0,
             'non-repeat': 0}
for line in total_file:
    for key in sv_rep_dict[line.split('\t')[4]]:
        total_rep[key] += 1

# Extract recall by each tool for each SV repeat family
def repeat_analysis(tp_file, sv_rep_dict, total_rep):
    tp_list = []
    with open(tp_file) as f:
        for line in f:
            svname = line.split('\t')[2].split('^')[0].split('-')[0].strip()
            tp_list.append(svname)
    tp_list = list(set(tp_list))
    data_rep = {"DNA": 0, "LINE": 0, "SINE": 0, "LTR": 0, "Low_complexity": 0, "Satellite": 0, "Simple_repeat": 0, 'others': 0,
                'non-repeat': 0}
    for i in tp_list:
        for key in sv_rep_dict[i]:
            data_rep[key] += 1
    for key in data_rep:
        data_rep[key] = round(data_rep[key]/total_rep[key], 3)
    return data_rep

# Execute repeat_analysis on each tool
nv1homo_rep = repeat_analysis(nv1homo, sv_rep_dict, total_rep)
snif1homo_rep = repeat_analysis(snif1homo, sv_rep_dict, total_rep)
sv1homo_rep = repeat_analysis(sv1homo, sv_rep_dict, total_rep)
ns1homo_rep = repeat_analysis(ns1homo, sv_rep_dict, total_rep)
pk1homo_rep = repeat_analysis(pk1homo, sv_rep_dict, total_rep)
nb1homo_rep = repeat_analysis(nb1homo, sv_rep_dict, total_rep)
de1homo_rep = repeat_analysis(de1homo, sv_rep_dict, total_rep)

nv1het_rep = repeat_analysis(nv1het, sv_rep_dict, total_rep)
snif1het_rep = repeat_analysis(snif1het, sv_rep_dict, total_rep)
sv1het_rep = repeat_analysis(sv1het, sv_rep_dict, total_rep)
ns1het_rep = repeat_analysis(ns1het, sv_rep_dict, total_rep)
pk1het_rep = repeat_analysis(pk1het, sv_rep_dict, total_rep)
nb1het_rep = repeat_analysis(nb1het, sv_rep_dict, total_rep)
de1het_rep = repeat_analysis(de1het, sv_rep_dict, total_rep)

# Plot barplots
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Patch

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
params = {'axes.labelsize': 9, 'axes.titlesize': 9, 'legend.fontsize': 7, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
matplotlib.rcParams.update(params)

df_homo = pd.DataFrame([nv1homo_rep, snif1homo_rep, sv1homo_rep, ns1homo_rep, pk1homo_rep, nb1homo_rep, de1homo_rep])
tools = ['nanovar', 'sniffles', 'svim', 'nanosv', 'picky', 'novobreak', 'delly']
df_homo["Tools"] = tools

df_het = pd.DataFrame([nv1het_rep, snif1het_rep, sv1het_rep, ns1het_rep, pk1het_rep, nb1het_rep, de1het_rep])
tools = ['nanovar', 'sniffles', 'svim', 'nanosv', 'picky', 'novobreak', 'delly']
df_het["Tools"] = tools

colors = ['mediumslateblue', '#23a7e4', 'saddlebrown', '#ffd919', '#3ca63b', 'crimson', 'darkorange']

# Create plot only for SINE and LINE repeat families
fig = plt.figure(figsize=(3.7, 3.9))
ax1 = fig.add_subplot(221)
sns.barplot(x="Tools", y="SINE", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('SINE')
ax12 = ax1.twinx()
ax1.set_ylim(0.0, 1.0)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax12.set_ylabel("No. of SVs", rotation=270)
ax12.set_ylim(0, total_rep['SINE'])
plt.yticks(rotation=45)

ax2 = fig.add_subplot(222)
sns.barplot(x="Tools", y="LINE", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('LINE')
ax22 = ax2.twinx()
ax2.set_ylim(0.0, 1.0)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax22.set_ylabel("No. of SVs", rotation=270)
ax22.set_ylim(0, total_rep['LINE'])
plt.yticks(rotation=45)

ax3 = fig.add_subplot(223)
sns.barplot(x="Tools", y="SINE", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('SINE')
ax32 = ax3.twinx()
ax3.set_ylim(0.0, 1.0)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax32.set_ylabel("No. of SVs", rotation=270)
ax32.set_ylim(0, total_rep['SINE'])
plt.yticks(rotation=45)

ax4 = fig.add_subplot(224)
sns.barplot(x="Tools", y="LINE", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('LINE')
ax42 = ax4.twinx()
ax4.set_ylim(0.0, 1.0)
ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax42.set_ylabel("No. of SVs", rotation=270)
ax42.set_ylim(0, total_rep['LINE'])
plt.yticks(rotation=45)

fig.tight_layout()
plt.savefig('fig.svg', format='svg', dpi=350)



# Create plot for the rest of repeat families (homozygous)
fig = plt.figure(figsize=(5, 3))
ax1 = fig.add_subplot(231)
sns.barplot(x="Tools", y="DNA", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('DNA transposon')
ax12 = ax1.twinx()
ax1.set_ylim(0.0, 1.0)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax12.set_ylabel("No. of SVs", rotation=270)
ax12.set_ylim(0, total_rep['DNA'])
plt.yticks(rotation=45)

ax2 = fig.add_subplot(232)
sns.barplot(x="Tools", y="LTR", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('LTR')
ax22 = ax2.twinx()
ax2.set_ylim(0.0, 1.0)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax22.set_ylabel("No. of SVs", rotation=270)
ax22.set_ylim(0, total_rep['LTR'])
plt.yticks(rotation=45)

ax3 = fig.add_subplot(233)
sns.barplot(x="Tools", y="Low_complexity", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Low complexity')
ax32 = ax3.twinx()
ax3.set_ylim(0.0, 1.0)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax32.set_ylabel("No. of SVs", rotation=270)
ax32.set_ylim(0, total_rep['Low_complexity'])
plt.yticks(rotation=45)

ax4 = fig.add_subplot(234)
sns.barplot(x="Tools", y="Satellite", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Satellite')
ax42 = ax4.twinx()
ax4.set_ylim(0.0, 1.0)
ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax42.set_ylabel("No. of SVs", rotation=270)
ax42.set_ylim(0, total_rep['Satellite'])
plt.yticks(rotation=45)

ax5 = fig.add_subplot(235)
sns.barplot(x="Tools", y="Simple_repeat", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Simple repeat')
ax52 = ax5.twinx()
ax5.set_ylim(0.0, 1.0)
ax5.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax52.set_ylabel("No. of SVs", rotation=270)
ax52.set_ylim(0, total_rep['Simple_repeat'])
plt.yticks(rotation=45)

ax6 = fig.add_subplot(236)
handle1 = sns.barplot(x="Tools", y="non-repeat", data=df_homo, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Non-repetitive')
ax62 = ax6.twinx()
ax6.set_ylim(0.0, 1.0)
ax6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax62.set_ylabel("No. of SVs", rotation=270)
ax62.set_ylim(0, total_rep['non-repeat'])
plt.yticks(rotation=45)

legend_elements = [Patch(facecolor='mediumslateblue', label='NanoVar'),
                   Patch(facecolor='#23a7e4', label='Sniffles'),
                   Patch(facecolor='saddlebrown', label='SVIM'),
                   Patch(facecolor='#ffd919', label='NanoSV'),
                   Patch(facecolor='#3ca63b', label='Picky'),
                   Patch(facecolor='crimson', label='novoBreak'),
                   Patch(facecolor='darkorange', label='Delly')]
legend1 = fig.legend(handles=legend_elements, bbox_to_anchor=(0.55, 0.35), ncol=3, fancybox=True)
legend1.get_frame().set_edgecolor('black')

fig.tight_layout()
plt.savefig('fig.svg', format='svg', dpi=350)



# Create plot for the rest of repeat families (heterozygous)
fig = plt.figure(figsize=(5, 3))
ax1 = fig.add_subplot(231)
sns.barplot(x="Tools", y="DNA", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('DNA transposon')
ax12 = ax1.twinx()
ax1.set_ylim(0.0, 1.0)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax12.set_ylabel("No. of SVs", rotation=270)
ax12.set_ylim(0, total_rep['DNA'])
plt.yticks(rotation=45)

ax2 = fig.add_subplot(232)
sns.barplot(x="Tools", y="LTR", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('LTR')
ax22 = ax2.twinx()
ax2.set_ylim(0.0, 1.0)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax22.set_ylabel("No. of SVs", rotation=270)
ax22.set_ylim(0, total_rep['LTR'])
plt.yticks(rotation=45)

ax3 = fig.add_subplot(233)
sns.barplot(x="Tools", y="Low_complexity", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Low complexity')
ax32 = ax3.twinx()
ax3.set_ylim(0.0, 1.0)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax32.set_ylabel("No. of SVs", rotation=270)
ax32.set_ylim(0, total_rep['Low_complexity'])
plt.yticks(rotation=45)

ax4 = fig.add_subplot(234)
sns.barplot(x="Tools", y="Satellite", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Satellite')
ax42 = ax4.twinx()
ax4.set_ylim(0.0, 1.0)
ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax42.set_ylabel("No. of SVs", rotation=270)
ax42.set_ylim(0, total_rep['Satellite'])
plt.yticks(rotation=45)

ax5 = fig.add_subplot(235)
sns.barplot(x="Tools", y="Simple_repeat", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Simple repeat')
ax52 = ax5.twinx()
ax5.set_ylim(0.0, 1.0)
ax5.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax52.set_ylabel("No. of SVs", rotation=270)
ax52.set_ylim(0, total_rep['Simple_repeat'])
plt.yticks(rotation=45)

ax6 = fig.add_subplot(236)
sns.barplot(x="Tools", y="non-repeat", data=df_het, palette=colors)
plt.xticks([], [])
plt.ylabel('Recall')
plt.title('Non-repetitive')
ax62 = ax6.twinx()
ax6.set_ylim(0.0, 1.0)
ax6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax62.set_ylabel("No. of SVs", rotation=270)
ax62.set_ylim(0, total_rep['non-repeat'])
plt.yticks(rotation=45)

fig.tight_layout()
plt.savefig('fig.svg', format='svg', dpi=350)

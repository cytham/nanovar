"""
Functions to calculate the upper limit for long-read SV depth coverage.

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pybedtools import BedTool
from scipy.interpolate import make_interp_spline, BSpline


# Generate upper overlap limit, depth of coverage, and coverage curve plot
def ovl_upper(total_gsize, contig_len_dict, basecov, subdata, wk_dir):
    # matplotlib.use('Agg')
    bed_str = normalbed(subdata)
    bed = BedTool(bed_str, from_string=True)
    bed = bed.sort()
    n = ngenerate(total_gsize)
    gsize_path = make_gsize(contig_len_dict, wk_dir)
    x = BedTool()
    ran_bed = x.random(l=100, n=n, seed=3, g=gsize_path)
    ran_bed = ran_bed.sort()
    insect = ran_bed.intersect(bed, wa=True, wb=True)
    cov_list = []
    for line in insect:
        line = list(line)
        cov_list.append(line[3])
    cov_dict = {}
    for key in cov_list:
        try:
            cov_dict[key] += 1
        except KeyError:
            cov_dict[key] = 1
    data = []
    for key in cov_dict:
        data.append(cov_dict[key])
    data = [float(i) for i in data]
    zerolist = [0.0]*(n-len(data))
    data2 = data + zerolist
    med = np.median(data2)
    medad = mad(data2)
    curve(data2, n, round((medad*6) + med, 0), wk_dir)
    depth = round(float(basecov)/total_gsize, 2)
    if depth < 1:
        maxovl = 10
    else:
        maxovl = round((medad * 4) + med, 1)
    return maxovl, depth


# Make genome size file
def make_gsize(contig_len_dict, wk_dir):
    path = os.path.join(wk_dir, 'genome.sizes')
    data = open(path, 'w')
    tmp = []
    for contig in contig_len_dict:
        checkcontigname(contig)
        tmp.append(contig + '\t' + str(contig_len_dict[contig]))
    data.write('\n'.join(tmp))
    data.close()
    return path


# Check contig name
def checkcontigname(contig):
    if "~" in contig or ":" in contig or "-" in contig:
        raise Exception("Contig name %s contains invalid symbols [ ~ : - ]" % contig)


# Generate number of genomic points
def ngenerate(gsize):
    if gsize > 10000:
        return 10000
    else:
        return int(gsize*0.9)


# Function to generate bed file from subdata
def normalbed(subdata):
    totalbed = []
    for line in subdata:
        coord1 = int(line.split('\t')[1])
        coord2 = int(line.split('\t')[1]) + int(line.split('\t')[2])
        totalbed.append(line.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + line.split('\t')[4])
    bed_str = '\n'.join(totalbed)
    return bed_str


# Median absolute deviation
def mad(x):
    b = 1.4826
    return b*np.mean(abs(x-np.median(x)))


# Plot curve
def curve(data, n, upper_limit, wk_dir):
    c = range(max(int(upper_limit), 10))
    p = []
    for i in c:
        p.append(data.count(i))
    p.append(n - sum(p))
    y = np.array([(float(z)/n) for z in p])
    # theoretical = [0.0915, 0.0441, 0.1032, 0.1498, 0.1739, 0.1626,
    # 0.1132, 0.0808, 0.0412, 0.0247, 0.0097, 0.0028, 0.0015,
    # 0.0006, 0.0002, 0.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    c = np.array(range(max(int(upper_limit), 10)+1))
    xnew = np.linspace(c.min(), c.max(), 100)
    spl = make_interp_spline(c, y)
    smooth = spl(xnew)
    params = {'axes.labelsize': 14, 'axes.titlesize': 17, 'legend.fontsize': 10, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
              'font.family': 'Arial, Helvetica, sans-serif'}
    matplotlib.rcParams.update(params)
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    ax.plot(xnew, smooth, color='#403f7d', linewidth=2.0)
    ax.set_facecolor('#ebebff')
    plt.text(int(upper_limit), y[-1], '>=' + str(int(upper_limit)))
    ax.grid(color='w', linestyle='-', linewidth=1)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.1%}'.format(x) for x in vals])
    plt.ylabel('Percentage')
    plt.xlabel('Depth of coverage')
    ax.legend(['Input sequencing data'], loc='upper right', ncol=1, fancybox=True)
    plt.savefig(os.path.join(wk_dir, 'fig', 'depth_of_coverage.png'),
                bbox_inches='tight',
                dpi=100,
                facecolor=fig.get_facecolor(),
                edgecolor='none')

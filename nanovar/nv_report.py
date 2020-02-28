"""
Functions to create HTML report.

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
import math
import datetime
import numpy as np
import nanovar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from distutils.dir_util import copy_tree
from nanovar import __version__


# Create HTML report
def create_report(wk_dir, contig_len_dict, thres, read_path, ref_path, rlen_dict, read_name, num_limit, ratio_limit):
    timenow = datetime.datetime.now()
    fwd = os.path.join(os.getcwd(), wk_dir, 'fig')
    fwd_fig = './fig'
    threshold = thres
    num_header = len(contig_len_dict) + 24
    vcf_path = os.path.join(wk_dir, '%s.nanovar.total.vcf' % read_name)
    vcf_path_pass = os.path.join(wk_dir, '%s.nanovar.pass.vcf' % read_name)
    vcf_data = open(vcf_path, 'r').read().splitlines()
    vcf = sorted(vcf_data[num_header:], key=lambda x: float(x.split('\t')[5]), reverse=True)
    # Creating variables
    scorelist, ratiolist, lcovlist = [], [], []
    svdict = {'DEL': 0, 'INV': 0, 'DUP': 0, 'INS': 0, 'BND': 0}
    dellen, inslen, invlen, duplen, bndlen = [], [], [], [], []
    delnolen, insnolen, invnolen, dupnolen, bndnolen = 0, 0, 0, 0, 0
    data = []
    totalsv, n = 0, 0
    for i in vcf:
        totalsv += 1
        scorelist.append(float(i.split('\t')[5]))
        lcov = int(i.split('\t')[9].split(':')[2].split(',')[1])
        lcovlist.append(lcov)
        ncov = int(i.split('\t')[9].split(':')[2].split(',')[0])
        ratio = round(lcov/(lcov + ncov), 3)
        ratiolist.append(ratio)
        if float(i.split('\t')[5]) >= threshold:
            if ratio <= ratio_limit:
                n += 1
                if i.split('\t')[4].startswith('<'):
                    sv_type = i.split('\t')[4].strip('<>')
                    sv_len = i.split('\t')[7].split(';')[2].split('=')[1]
                    data.append([str(n), sv_type, i.split('\t')[0], i.split('\t')[1],
                                 i.split('\t')[4], i.split('\t')[7].split(';')[1].split('=')[1],
                                 sv_len, i.split('\t')[5],
                                 str(lcov), str(ncov), str(ratio), i.split('\t')[9].split(':')[0],
                                 i.split('\t')[2]])
                    svdict[sv_type] += 1
                    if sv_type == 'DEL':
                        try:
                            dellen.append(svlencap(int(sv_len)))
                        except ValueError:
                            delnolen += 1
                    elif sv_type == 'INS':
                        try:
                            inslen.append(svlencap(int(sv_len)))
                        except ValueError:
                            insnolen += 1
                    elif sv_type == 'INV':
                        try:
                            invlen.append(svlencap(int(sv_len)))
                        except ValueError:
                            invnolen += 1
                    elif sv_type == 'DUP':
                        try:
                            duplen.append(svlencap(int(sv_len)))
                        except ValueError:
                            dupnolen += 1
                else:
                    sv_type = i.split('\t')[7].split(';')[5].split('=')[1]
                    chrom2 = i.split('\t')[4].split(':')[0].strip('N[]')
                    pos2 = i.split('\t')[4].split(':')[1].strip('N[]')
                    sv_len = i.split('\t')[7].split(';')[2].split('=')[1]
                    data.append([str(n), sv_type, i.split('\t')[0], i.split('\t')[1],
                                 chrom2, pos2, sv_len, i.split('\t')[5], str(lcov), str(ncov),
                                 str(ratio), i.split('\t')[9].split(':')[0], i.split('\t')[2]])
                    svdict['BND'] += 1
                    try:
                        bndlen.append(svlencap(int(sv_len)))
                    except ValueError:
                        bndnolen += 1
    totalsvlen = [dellen, inslen, invlen, bndlen, duplen]
    tab = [[len(dellen), len(inslen), len(invlen), len(bndlen), len(duplen)], [delnolen, insnolen, invnolen, bndnolen, dupnolen]]
    # Setting global figure parameters
    params = {'axes.labelsize': 14, 'axes.titlesize': 17, 'legend.fontsize': 10,
              'xtick.labelsize': 12, 'ytick.labelsize': 12, 'font.family': 'Arial, Helvetica, sans-serif'}
    matplotlib.rcParams.update(params)
    # Make plots
    scatter_plots(fwd, scorelist, ratiolist, lcovlist, threshold)
    sv_len_dict(fwd, totalsvlen, tab)
    sv_type_dist(svdict, fwd)
    read_len_dict(rlen_dict, fwd)
    # Resize data for table
    data2 = data[0:num_limit]
    # Write HTML
    create_html(data2, fwd_fig, wk_dir, vcf_path_pass, timenow, read_name, read_path, ref_path, threshold, n, totalsv)
    # Copy css and js directories
    css = os.path.join(os.path.dirname(nanovar.__file__), 'css')
    js = os.path.join(os.path.dirname(nanovar.__file__), 'js')
    css_to = os.path.join(os.getcwd(), wk_dir, 'css')
    js_to = os.path.join(os.getcwd(), wk_dir, 'js')
    null = copy_tree(css, css_to)
    null = copy_tree(js, js_to)


# SV Length cap
def svlencap(x):
    if x >= 10000:
        return 10000
    else:
        return x


# Scatter plots generation
def scatter_plots(fwd, scorelist, ratiolist, lcovlist, threshold):
    # Scatter plot Score vs Read ratio
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    ax.scatter(ratiolist, scorelist, c='#7d3f5d', alpha=0.1)
    ax.axhline(y=threshold, linewidth=1, color='firebrick')
    ax.annotate("Threshold=" + str(threshold), xy=(0.8, threshold+0.2))
    ax.set_facecolor('#ebebff')
    plt.ylabel('Confidence score')
    plt.xlabel('Breakend read ratio')
    plt.ylim(bottom=-0.3)
    plt.savefig(os.path.join(fwd, 'scatter1.png'), bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')
    # Scatter plot Score vs Supporting reads
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    if lcovlist:
        mcov = max(max(lcovlist), 10)
    else:
        mcov = 10
    ax.scatter(lcovlist, scorelist, c='#3f5d7d', alpha=0.1)
    ax.axhline(y=threshold, linewidth=1, color='firebrick')
    ax.annotate("Threshold=" + str(threshold), xy=(mcov - 2.8, threshold+0.2))
    ax.set_facecolor('#ebebff')
    plt.ylabel('Confidence score')
    plt.xlabel('Number of breakend-supporting reads')
    plt.ylim(bottom=-0.3)
    plt.xlim(right=mcov)
    plt.savefig(os.path.join(fwd, 'scatter2.png'), bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')


# Boxplot SV length distribution
def sv_len_dict(fwd, totalsvlen, tab):
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    label = ["Deletion", "Insertion", "Inversion", "Breakend", "TandemDup"]
    colors = ['#7d3f5d', '#3f5d7d', '#7d5f3f', '#403f7d', '#3f7c7d']
    medianprop = dict(linewidth=1.5, color='black')
    bp = ax.boxplot(totalsvlen, whis=[5, 95], medianprops=medianprop, showfliers=False, patch_artist=True)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    ax.table(cellText=tab, rowLabels=['Known size', 'Unknown size'], colLabels=label, loc='bottom')
    ax.set_facecolor('#ebebff')
    ax.xaxis.grid(False)
    ax.yaxis.grid(color='white', linewidth=1)
    plt.ylabel('Length (base pair)')
    plt.savefig(os.path.join(fwd, 'sv_lengths.png'), bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(),
                edgecolor='none')


# Donut SV types
def sv_type_dist(svdict, fwd):
    fig, ax = plt.subplots(figsize=(10, 5), subplot_kw=dict(aspect="equal"))
    fig.patch.set_facecolor('#f6f7f9')
    svnamedict = {"DEL": " Deletions", "INS": " Insertions",
                  "INV": " Inversions", "DUP": " TandemDups",
                  "BND": " Breakends (TLO/TPO)"}
    label = []
    data = []
    for key in svdict:
        if svdict[key] > 0:
            label.append(str(svdict[key]) + svnamedict[key])
            data.append(svdict[key])

    def func(pct, allvals):
        absolute = int(pct/100.*np.sum(allvals))
        return "{:.1f}%".format(pct, absolute)

    explode = []
    for i in range(len(data)):
        explode.append(0.02)
    explode = tuple(explode)  # explode = (0.02, 0.02, 0.02, 0.02, 0.02)
    wedges, text, autotexts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=-40, autopct=lambda pct: func(pct, data),
                                     textprops=dict(color="white"), counterclock=False, pctdistance=0.75,
                                     colors=['#7d3f5d', '#7d5f3f', '#3f5d7d', '#3f7c7d', '#403f7d'],
                                     explode=explode)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        if data[i] > 0:
            ax.annotate(label[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                        horizontalalignment=horizontalalignment, **kw)

    plt.setp(autotexts, size=11)
    ax.text(0, 0, 'SV types', ha='center')
    plt.savefig(os.path.join(fwd, 'sv_type_donut.png'), bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(),
                edgecolor='none')


# Create read length distribution plot
def read_len_dict(rlen_dict, fwd):
    qlen = measureqlen(rlen_dict)
    m = math.ceil(max(np.log10(np.array(qlen))))
    e = math.floor(min(np.log10(np.array(qlen))))
    bins = 10**np.linspace(e, m, m*20)
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    ax.hist(qlen, bins=bins, color="#403f7d", edgecolor='black', linewidth=0.5)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.ticklabel_format(useOffset=False, style='plain')
    ax.set_facecolor('#ebebff')
    plt.xlabel("Read length (bases)")
    plt.ylabel("Number of reads")
    plt.savefig(os.path.join(fwd, 'read_length_dist.png'), bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(),
                edgecolor='none')


# Get read length distribution
def measureqlen(rlen_dict):
    qlen = []
    for key in rlen_dict:
        qlen.append(rlen_dict[key])
    return qlen


# Create html
def create_html(data, fwd, wk_dir, vcf_path, timenow, read_name, read_path, ref_path, threshold, n, totalsv):
    html = open(os.path.join(wk_dir, '%s.nanovar.pass.report.html' % read_name), 'w')
    begin = """<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <meta http-equiv="x-ua-compatible" content="ie=edge">
        <title>NanoVar Report</title>
        <!-- Font Awesome -->
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
        <!-- Bootstrap core CSS -->
        <link href="css/bootstrap.min.css" rel="stylesheet">
        <!-- Material Design Bootstrap -->
        <link href="css/mdb.min.css" rel="stylesheet">
        <!-- DataTable CSS -->
        <link href="css/jquery.dataTables.min.css" rel="stylesheet">
        <!-- DataTable buttons CSS  -->
        <link href="css/buttons.dataTables.min.css" rel="stylesheet">
        <style>
            body {
                background-color: #f6f7f9;
            }
            h1 { 
                margin-left: 4px;
                font-weight: bold;
                font-size: 40px;
                font-family: Arial, Helvetica, sans-serif;
            }
            h2, h3, h4, h5, p { 
                margin-left: 4px;
                font-family: Arial, Helvetica, sans-serif;
            }
            table {
                font-family: Arial, Helvetica, sans-serif;
            }
            #image {
                text-align:center;
            }
            #margin {
                margin-left: 35px;
                margin-right: 35px;
                margin-top: 25px;
            }
            .floatleft {
                float: left;
                width: 45%;
                margin-right: 5%;
                padding-top: 80px
            }
            .floatright {
                float: right;
                width: 45%;
                margin-left: 5%;
            }
            .container {
                overflow;
                max-width: 1400px;
            }
            .aligncenter {
                float: center;
            }
            .alignleft {
                float: left;
            }
            .alignright {
                float: right;
            }
        </style>
    </head>
    <body>
    <div id="margin">
        <div>
            <h1 style="text-align:center;">NanoVar Report</h1>
            <h5 style="text-align:center;">Version: NanoVar """ + __version__ + """</h5>
            <h5 style="text-align:center;">""" + str(timenow.strftime("%a %d %B %Y")) + """</h5>
            <h5 style="text-align:center;">""" + read_name + """</h5>
            <br>
            <br>
            <div style="background-color:#3f5d7d;color:#f6f7f9">
                <h4 style="text-align:center;">Run details</h4>
            </div>
        </div>
        <div style="clear: both;"></div>
        <br>
        <table width="50%" class="table-sm table-bordered" style="table-layout:fixed" align="center">
            <col width="150">
            <thead>
            <tr>
                <th style="background-color:#2e6096; color:#f6f7f9; font-weight: bold; text-align:center;">Variable</th>
                <th style="background-color:#2e6096; color:#f6f7f9; font-weight: bold; text-align:center;">Value</th>
            </tr>
            </thead>
            <tbody>
            <tr>
                <td>Input reads:</td>
                <td style="word-wrap:break-word;">""" + read_path + """</td>
            </tr>
            <tr>
                <td>Reference genome:</td>
                <td style="word-wrap:break-word;">""" + ref_path + """</td>
            </tr>
            <tr>
                <td>Output VCF:</td>
                <td style="word-wrap:break-word;">""" + vcf_path + """</td>
            </tr>
            <tr>
                <td>SV score threshold:</td>
                <td>""" + str(threshold) + """</td>
            </tr>
            </tbody>
        </table>
        <br>
        <br>
        <br>
        <br>
        <div style="background-color:#366b6c;color:#f6f7f9;">
            <h3 style="text-align:center;">Results</h3>
        </div>
        <p style="clear: both;">
        <br>
        <h4 style="text-align:center;"><u>1. Table of output SVs</u></h4>
        <h5 style="text-align:center;">Showing """ + str(len(data)) + """ out of """ + str(totalsv) + """ total SVs</h5>
        <table id="NanoVar_report_table" class="table table-striped table-bordered table-sm" cellspacing="0" width="100%">
            <thead>
                <tr>
                    <th class="th-sm">#
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">SV type
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Chrom1
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Pos1
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Chrom2
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Pos2
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Length
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Score
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">No. of breakend-supporting reads
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">No. of breakend-opposing reads
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th title="No. of breakend-supporting reads/total reads at breakend" class="th-sm">Breakend read ratio
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th class="th-sm">Approximated genotype
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                    <th title="Last five digits define the breakend ID code" class="th-sm">SV ID
                        <i class="fa fa-sort float-right" aria-hidden="true"></i>
                    </th>
                </tr>
            </thead>
            <tbody>
    """
    html.write(begin)

    # Create table body
    row = ""
    for i in data:
        row += """          <tr>
    """
        for j in i:
            row += "                <td>" + j + """</td>
    """
        row += """          </tr>
    """
    row += """      </tbody>
        </table>
        <br>
        <br>
        <div id="image">
            <figure>
                <h4 style="text-align:center;"><u>2. Distribution of SV types</u></h4>
                <img src=""" + '"' + fwd + """/sv_type_donut.png" alt="2. Distribution of SV types" title=""" + '"' + fwd + \
           '/sv_type_donut.png' + '"' + """>
            </figure>
            <br>
            <br>
            <figure>
                <h4 style="text-align:center;"><u>3. Size distribution of SVs across SV types</u></h4>
                <img src=""" + '"' + fwd + """/sv_lengths.png" alt="3. Size distribution of SVs across SV types" title=""" + \
           '"' + fwd + '/sv_lengths.png' + '"' + """>
            </figure>
            <br>
            <br>
            <figure>
                <h4 style="text-align:center;"><u>4. Scatter plot between SV confidence score and read ratio</u></h4>
                <img src=""" + '"' + fwd + """/scatter1.png" alt="4. Scatter plot between SV confidence score and read ratio" 
                title=""" + '"' + fwd + '/scatter1.png' + '"' + """>
            </figure>
            <br>
            <br>
            <figure>
                <h4 style="text-align:center;"><u>5. Scatter plot between SV confidence score and read depth</u></h4>
                <img src=""" + '"' + fwd + """/scatter2.png" alt="5. Scatter plot between SV confidence score and read depth" 
                title=""" + '"' + fwd + '/scatter2.png' + '"' + """>
            </figure>
            <br>
            <br>
            <div style="background-color:#403f7d;color:#f6f7f9;">
                <h3 style="text-align:center;">QC Statistics</h3>
            </div>
            <br>
            <figure>
                <h4 style="text-align:center;"><u>6. Distribution of read lengths of input FASTA</u></h4>
                <img src=""" + '"' + fwd + """/read_length_dist.png" alt="6. Distribution of read lengths of input FASTA" 
                title=""" + '"' + fwd + '/read_length_dist.png' + '"' + """>
            </figure>
            <br>
            <br>
            <figure>
                <h4 style="text-align:center;"><u>7. Depth of coverage across genome</u></h4>
                <img src=""" + '"' + fwd + """/depth_of_coverage.png" alt="7. Depth of coverage across genome" title=""" + '"' \
           + fwd + '/figures/depth_of_coverage.png' + '"' + """>
            </figure>
            <br>
        </div>
        <br>
        <br>
        <!-- SCRIPTS -->
        <!-- JQuery -->
        <script type="text/javascript" src="js/jquery-3.3.1.min.js"></script>
        <!-- datatable javascript  -->
        <script type="text/javascript" src="js/jquery.dataTables.min.js"></script>
        <!-- datatable buttons javascript  -->
        <script type="text/javascript" src="js/dataTables.buttons.min.js"></script>
        <!-- buttons flash javascript  -->
        <script type="text/javascript" src="js/buttons.flash.min.js"></script>
        <!-- buttons jszip javascript  -->
        <script type="text/javascript" src="js/jszip.min.js"></script>
        <!-- buttons html5 buttons javascript  -->
        <script type="text/javascript" src="js/buttons.html5.min.js"></script>
        <!-- Bootstrap tooltips -->
        <script type="text/javascript" src="js/popper.min.js"></script>
        <!-- Bootstrap core JavaScript -->
        <script type="text/javascript" src="js/bootstrap.min.js"></script>
        <!-- MDB core JavaScript -->
        <script type="text/javascript" src="js/mdb.min.js"></script>
        <script type="text/javascript">
            $('#NanoVar_report_table').DataTable({
                "scrollX": true,
                "scrollY": 200,
                dom: 'Bfrtip',
                buttons: [
                    'copy', 'csv', 'excel'
                ]
            });
            $('.dataTables_length').addClass('bs-select');
        </script>
    </div>
    </body>
    </html>
    """
    html.write(row)
    html.close()

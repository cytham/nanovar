"""
Functions for verifying files.

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


import gzip
import logging


def fastx_valid(reads, ext):
    c = 0
    lineskip = 0
    totalc = 0
    corrupted = 0
    # rlen_dict = {}
    if ext == "gz":
        with gzip.open(reads, 'rt') as f:
            line1n2 = [next(f) for x in range(2)]
            try:
                line3n4 = [next(f) for x in range(2)]
            except StopIteration:
                pass
            if "line3n4" in locals():
                if line1n2[0][0] == ">" and line3n4[0][0] == ">":
                    xformat = "fasta"
                    # rlen_dict[line1n2[0].split()[0][1:].strip()] = len(line1n2[1].strip())
                    # rlen_dict[line3n4[0].split()[0][1:].strip()] = len(line3n4[1].strip())
                    for line in f:
                        c += 1
                        if lineskip == 0:
                            if line[0] != ">":
                                corrupted = 1
                                break
                            rname = line.split()[0][1:].strip()
                            lineskip = 1
                        else:
                            # rlen_dict[rname] = len(line.strip())
                            lineskip = 0
                            if c == 100000:
                                totalc = totalc + c
                                c = 0
                elif line1n2[0][0] == "@" and line3n4[0][0] == "+":
                    xformat = "fastq"
                    # rlen_dict[line1n2[0].split()[0][1:].strip()] = len(line1n2[1].strip())
                    for line in f:
                        c += 1
                        if lineskip == 0:
                            if line[0] != "@":
                                corrupted = 1
                                break
                            rname = line.split()[0][1:].strip()
                            lineskip = 1
                        elif lineskip == 1:
                            # rlen_dict[rname] = len(line.strip())
                            lineskip = 2
                        elif lineskip == 2:
                            if line[0] != "+":
                                corrupted = 1
                                break
                            lineskip = 3
                        else:
                            lineskip += 1
                            if lineskip == 4:
                                lineskip = 0
                            if c == 100000:
                                totalc = totalc + c
                                c = 0
                else:
                    c = 1
                    corrupted = 1
            else:
                if line1n2[0][0] == ">":
                    xformat = "fasta"
                    # rlen_dict[line1n2[0].split()[0][1:].strip()] = len(line1n2[1].strip())
                    c = 0
                else:
                    c = 1
                    corrupted = 1
    elif ext == "txt":
        with open(reads) as f:
            line1n2 = [next(f) for x in range(2)]
            try:
                line3n4 = [next(f) for x in range(2)]
            except StopIteration:
                pass
            if "line3n4" in locals():
                if line1n2[0][0] == ">" and line3n4[0][0] == ">":
                    xformat = "fasta"
                    # rlen_dict[line1n2[0].split()[0][1:].strip()] = len(line1n2[1].strip())
                    # rlen_dict[line3n4[0].split()[0][1:].strip()] = len(line3n4[1].strip())
                    for line in f:
                        c += 1
                        if lineskip == 0:
                            if line[0] != ">":
                                corrupted = 1
                                break
                            rname = line.split()[0][1:].strip()
                            lineskip = 1
                        else:
                            # rlen_dict[rname] = len(line.strip())
                            lineskip = 0
                            if c == 100000:
                                totalc = totalc + c
                                c = 0
                elif line1n2[0][0] == "@" and line3n4[0][0] == "+":
                    xformat = "fastq"
                    # rlen_dict[line1n2[0].split()[0][1:].strip()] = len(line1n2[1].strip())
                    for line in f:
                        c += 1
                        if lineskip == 0:
                            if line[0] != "@":
                                corrupted = 1
                                break
                            rname = line.split()[0][1:].strip()
                            lineskip = 1
                        elif lineskip == 1:
                            # rlen_dict[rname] = len(line.strip())
                            lineskip = 2
                        elif lineskip == 2:
                            if line[0] != "+":
                                corrupted = 1
                                break
                            lineskip = 3
                        else:
                            lineskip += 1
                            if lineskip == 4:
                                lineskip = 0
                            if c == 100000:
                                totalc = totalc + c
                                c = 0
                else:
                    c = 1
                    corrupted = 1
            else:
                if line1n2[0][0] == ">":
                    xformat = "fasta"
                    # rlen_dict[line1n2[0].split()[0][1:].strip()] = len(line1n2[1].strip())
                    c = 0
                else:
                    c = 1
                    corrupted = 1
    # Summarize and count reads
    if corrupted == 0:
        if xformat == 'fasta':
            if c == 0:
                return ["Pass", int(((totalc + c) / 2) + 1)]  # [Status, read count]
            elif c % 2 == 0:
                return ["Pass", int(((totalc + c) / 2) + 2)]  # [Status, read count]
            else:
                return ["Fail", int(totalc + c) + 2]  # [Status, line error]
        elif xformat == 'fastq':
            if c % 4 == 0:
                return ["Pass", int(((totalc + c) / 4) + 1)]  # [Status, read count]
            else:
                return ["Fail", int(totalc + c) + 4]  # [Status, line error]
    elif corrupted == 1:
        return ["Fail", int(totalc + c) + 4]  # [Status, line error]


def bed_valid(filter_path, contig_len_dict):

    # Make contig BED dict
    contig_bed_dict = {}
    contig_list = []
    for cid in contig_len_dict:
        contig_bed_dict[cid] = [0, contig_len_dict[cid]]
        contig_list.append(cid)

    # Validate BED file to genome
    c = 0
    with open(filter_path) as bed:
        for line in bed:
            c += 1
            cid = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            if cid in contig_list:
                try:
                    start = int(start)
                    end = int(end)
                except TypeError:
                    logging.critical("Error: Genome filter BED non-integer at line %s" % str(c))
                    raise TypeError("Error: Genome filter BED non-integer at line %s" % str(c))
                if contig_bed_dict[cid][0] <= start and contig_bed_dict[cid][1] >= end:
                    if start <= end:
                        return True
                    else:
                        logging.critical("Error: Genome filter BED start>end at line %s" % str(c))
                        raise TypeError("Error: Genome filter BED start>end at line %s" % str(c))
                else:
                    logging.critical("Error: Genome filter BED outside contig boundaries at line %s" % str(c))
                    raise TypeError("Error: Genome filter BED outside contig boundaries at line %s" % str(c))
            else:
                logging.critical("Error: Genome filter BED contig id not found in reference at line %s" % str(c))
                raise TypeError("Error: Genome filter BED contig id not found in reference at line %s" % str(c))

"""
Functions to verify input BED files.

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


import logging


def bed_valid(filter_path, contig_len_dict):

    # Make contig BED dict
    contig_bed_dict = {}
    contig_list = []
    for id in contig_len_dict:
        contig_bed_dict[id] = [0, contig_len_dict[id]]
        contig_list.append(id)

    # Validate BED file to genome
    c = 0
    with open(filter_path) as bed:
        for line in bed:
            c += 1
            id = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            if id in contig_list:
                try:
                    start = int(start)
                    end = int(end)
                except TypeError:
                    logging.critical("Error: Genome filter BED non-integer at line %s" % str(c))
                    raise TypeError("Error: Genome filter BED non-integer at line %s" % str(c))
                if contig_bed_dict[id][0] < start and contig_bed_dict[id][1] > end:
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

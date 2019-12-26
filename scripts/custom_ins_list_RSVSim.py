import math
import numpy as np
import sys
from sys import argv
import random

if len(argv) != 7:
    sys.exit("Usage: python custom_ins_list_RSVSim.py viral.fa insertions.csv seed_no no_of_viral_ins_per_set no_of_sets output_path")

# Define paths
viral = argv[1]
hg = argv[2]
seed = int(argv[3])
# Set number of novel insertions
no_ins = int(argv[4])
no_iter = int(argv[5])
out_dir = argv[6]
np.random.seed(seed)


class InsertViral:
    
    def __init__(self, viral_fasta, ins_regions, nov, itera, out):
        self.nov = nov
        self.itera = itera
        self.out = out
        self.viralnames = []
        self.viralsize_dict = {}
        self.viralseq_dict = {}
        self.viralcounter = {}
        with open(viral_fasta) as data:
            for line in data:
                self.viralnames.append(line.split(' ')[0])
                fasta = next(data).strip()
                self.viralseq_dict[line.split(' ')[0]] = str(fasta)
                self.viralsize_dict[line.split(' ')[0]] = len(fasta)
                self.viralcounter[line.split(' ')[0]] = 1
        self.viralnum = len(self.viralnames)
        data = open(ins_regions, 'r').read().splitlines()
        # del data[0]
        self.coorlist = []
        for line in data:
            if line.split('\t')[0] != "Name":
                cline = line.split('\t')
                self.coorlist.append([cline[1], cline[2], cline[3], cline[4], cline[5], cline[6], cline[7]])
        random.shuffle(self.coorlist)

    def insert(self):
        # Check enough coordinates
        if len(self.coorlist) < (self.nov * self.itera):
            raise IndexError('Provided insertions.csv file insufficient for number of sets.')
        last = 0
        range_chunks = []
        for n in range(self.itera):
            range_chunks.append(range(last, last+self.nov))
            last += self.nov
        for n in range(self.itera):
            viralcount = self.viralcounter.copy()
            out = open(self.out + '/viral_ins_regions' + str(n+1) + '.csv', 'w')
            for i in range_chunks[n]:
                line = self.coorlist[i]
                name, start, end = self.viral_select(int(line[6]), viralcount)
                out.write(str(name) + '\t' + str(start) + '\t' + str(end) + '\t' + str(line[6]) + '\t' + str(line[3]) + '\t' +
                          str(line[4]) + '\t' + str(line[5]) + '\t' + str(line[6]) + '\n')
            out.close()
            out_vfasta = open(self.out + '/viral_ins_fasta' + str(n + 1) + '.fa', 'w')
            for i in self.viralnames:
                out_vfasta.write(i + "\n")
                dup_required = math.ceil(viralcount[i]/self.viralsize_dict[i]) + 1
                out_vfasta.write(str(self.viralseq_dict[i]*dup_required) + "\n")
            out_vfasta.close()

    def viral_select(self, size, viralcount):
        vindex = np.random.randint(0, self.viralnum)
        name = self.viralnames[vindex]
        start = viralcount[name]
        end = start + size - 1
        viralcount[name] += size + 10
        return name.strip('>'), start, end


if __name__ == '__main__':
    workspace = InsertViral(viral, hg, no_ins, no_iter, out_dir)
    workspace.insert()

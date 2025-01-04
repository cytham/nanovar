"""
Generate SV-supporting reads BAM file. 

@author: asmaa, cytham

# Test
wk_dir = '.'
vcf_file = 'sample.nanovar.pass.vcf'
sv_supp_tsv_file = 'sv_support_reads.tsv'
in_bam_path = 'input.bam'
input_type = 'bam'
ref_path = 'ref.fa'
vcf = os.path.join(wk_dir, vcf_file)
sv_supp_tsv = os.path.join(wk_dir, sv_supp_tsv_file)
create_sv_supp_bam(vcf, sv_supp_tsv, in_bam_path, wk_dir, input_type, ref_path)

"""

import os
import pandas as pd
import pysam

def create_sv_supp_bam(vcf, sv_supp_tsv, in_bam_path, wk_dir, input_type, ref_path):
    pass_ids = parse_pass_sv(vcf)
    supp_dict = parse_supp_tsv(sv_supp_tsv, pass_ids)
    in_bam = input_type_func(in_bam_path, input_type, ref_path)
    out_bam = os.path.join(wk_dir, "sv_support_reads.bam")
    with pysam.AlignmentFile(out_bam, "wb", template=in_bam) as bam_out:
        for seg in in_bam:
            try:
                id = supp_dict[seg.query_name]  # If segment is SV-supporting read
                seg.set_tag('nv', ','.join(id), value_type='Z', replace=True)  # Set SV-ID tag
                bam_out.write(seg)
            except KeyError:
                pass
    bam_out.close()

def parse_supp_tsv(sv_sup_tsv, pass_ids):
    supp_df = pd.read_csv(sv_sup_tsv, sep='\t', header=0)
    supp_df.rename(columns = {supp_df.columns[1]: 'readID'}, inplace = True)
    supp_dict = {}
    for index, row in supp_df.iterrows():
        for _id in row['readID'].split(','):
            try:
                _ = pass_ids[row['SV-ID']]  # If passed SV-ID
                id = _id.split('~')[0]
                try:
                    supp_dict[id].add(row['SV-ID'])
                except KeyError:
                    supp_dict[id] = set()
                    supp_dict[id].add(row['SV-ID'])
            except KeyError:
                pass
    return supp_dict

def parse_pass_sv(vcf):
    pass_ids = {}
    with open(vcf) as v:
        for line in v:
            if not line.startswith('#'):
                pass_ids[line.split('\t')[2]] = 0
    return pass_ids

def input_type_func(file_path, input_type, ref_path):
    if input_type in ['bam', 'raw']:
        sam = pysam.AlignmentFile(file_path, "rb")
    elif input_type == 'cram':
        sam = pysam.AlignmentFile(file_path, "rc", reference_filename=ref_path)
    return sam

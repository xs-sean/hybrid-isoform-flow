#!/usr/bin/python

import sys
import os
from decimal import Inexact

### Generate output gpd file with fusion pairs that have exon match with an annotated isoform 
### outfile.all will include all annotated matching alignments
### Requirements: pairs should be consecutive in input file

if len(sys.argv) >= 3:
    fusion_gpd_filename = sys.argv[1]
    ref_gpd_filename = sys.argv[2]
    output_gpd_filename = sys.argv[3]
else:
    sys.exit(1)


### Parses reference gpd file to generate list of transcripts sequence junction
### Updates transseq_dict, transname_dict dictionaries
def parse_ref_gpd_file(ref_gpd_filename, transseq_dict, transname_dict):
    
    ref_gpd_file = open(ref_gpd_filename, 'r')
    for line in ref_gpd_file:
        fields = line.split('\t')
        chrname = fields[2]
        gname = fields[0]
        transname = fields[1]
        strand = fields[3]
        startpos = fields[9].split(',')[1:-1]
        endpos = fields[10].split(',')[:-2]
        
        if (not transseq_dict.has_key(chrname)):
            transseq_dict[chrname] = []
            transname_dict[chrname] = dict()
        
        juncseq = []
        for idx in range(len(startpos)):
            juncseq.append(endpos[idx])
            juncseq.append(startpos[idx])
            
            juncseq_str = "-".join(juncseq)
            transseq_dict[chrname].append(juncseq_str)
            if (not transname_dict[chrname].has_key(juncseq_str)):
                transname_dict[chrname][juncseq_str] = []
            transname_dict[chrname][juncseq_str].append(":".join([gname, transname, strand]))
                    
        juncseq = []
        for idx in range(len(startpos)-1, -1, -1):
            juncseq.append(startpos[idx])
            juncseq.append(endpos[idx])
            
            juncseq_str = "-".join(list(reversed(juncseq)))
            transseq_dict[chrname].append(juncseq_str)
            if (not transname_dict[chrname].has_key(juncseq_str)):
                transname_dict[chrname][juncseq_str] = []
            transname_dict[chrname][juncseq_str].append(":".join([gname, transname, strand]))
    ref_gpd_file.close()


################################################################################

transseq_dict = dict()
transname_dict = dict()
parse_ref_gpd_file(ref_gpd_filename, transseq_dict, transname_dict)
juncseq_list = dict()
for chrname in transseq_dict.keys():
    juncseq_list[chrname] = list()
    for seq in transseq_dict[chrname]:
        juncseq_list[chrname].append(seq)
        
fusion_gpd_file = open(fusion_gpd_filename,'r')
output_gpd_file = open(output_gpd_filename,'w')
output_gpd_all_file = open(output_gpd_filename + '.all','w')
pair_idx = 0

while (True):
    line = fusion_gpd_file.readline()
    if not line:
        break
    if (pair_idx == 0):
        line_str = ['', '']
    
    fields = line.strip().split('\t')
    
    chrname = fields[2]
    gname = fields[0]
    transname = fields[1]
        
    startpos = fields[9].split(',')[1:-1]
    endpos = fields[10].split(',')[:-2]
    
    juncseq = []
    for idx in range(len(startpos)):
        juncseq.append(endpos[idx])
        juncseq.append(startpos[idx])
            
    juncseq_str = "-".join(juncseq)
    if (juncseq_str in juncseq_list[chrname]):
        fields[1] = "|".join(transname_dict[chrname][juncseq_str])
        fields[3] = transname_dict[chrname][juncseq_str][0].split(':')[-1]
        line_str[pair_idx] = "\t".join(fields + [juncseq_str]) + "\n"
        output_gpd_all_file.write(line_str[pair_idx])
    
    if (pair_idx == 1):
        if (line_str[0] != '' and line_str[1] != ''):
            output_gpd_file.write(line_str[0])
            output_gpd_file.write(line_str[1])
            
    
    pair_idx = 1- pair_idx
        
        
output_gpd_file.close()
fusion_gpd_file.close()
output_gpd_all_file.close()

            
    
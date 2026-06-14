#!/usr/bin/python

import sys
import os
import bisect
import re

gpd_filename = sys.argv[1]
output_filename = sys.argv[2]
ref_filename = sys.argv[3]

###
##############
def parse_ref_file(ref_filename):
    
    ref_file = open(ref_filename, 'r')
    gname_list = dict()
    gene_start_pos_list = dict()
    gene_end_pos_list = dict()
    gene_max_end_pos_list = dict()
        
    for line in ref_file:
        fields = line.split()
        gname = fields[0]
        start_pos = int(fields[4])
        end_pos = int(fields[5])
        rname = fields[2]
        idx = 0
        update_flag = True = False
        if (not gname_list.has_key(rname)):
            gname_list[rname] = []
            gene_start_pos_list[rname] = []
            gene_end_pos_list[rname] = []
        while (gname in gname_list[rname][idx:]):
            gname_idx =  gname_list[rname][idx:].index(gname) + idx
            idx = gname_idx + 1
            # Find an overlap
            if (min(gene_end_pos_list[rname][gname_idx], end_pos)  -
                max(gene_start_pos_list[rname][gname_idx], start_pos) > 0):
                gene_end_pos_list[rname][gname_idx] = max(gene_end_pos_list[rname][gname_idx], end_pos)
                gene_start_pos_list[rname][gname_idx] = min(gene_start_pos_list[rname][gname_idx], start_pos)
                update_flag = True
                break
        if not update_flag:
            gene_end_pos_list[rname].append(end_pos)
            gene_start_pos_list[rname].append(start_pos)
            gname_list[rname].append(gname)
            
    for chr in gname_list.keys():
        index_list = sorted(range(len(gene_start_pos_list[chr])),key=lambda x:gene_start_pos_list[chr][x])
        gname_list[chr] = [gname_list[chr][i] for i in index_list]
        gene_start_pos_list[chr] = [gene_start_pos_list[chr][i] for i in index_list]
        gene_end_pos_list[chr] = [gene_end_pos_list[chr][i] for i in index_list]
        gene_max_end_pos_list[chr] = [max(gene_end_pos_list[chr][:i]) for i in range(1, len(index_list)+1)]
    return [gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list]
            
###
################
def get_gnames(rname, tstart, tend,
               gname_list, gene_start_pos_list, 
               gene_end_pos_list, gene_max_end_pos_list):
    
    gnames = set()
    if (not gene_start_pos_list.has_key(rname)):
        return gnames
    index = min(bisect.bisect_left(gene_start_pos_list[rname], tend), len(gene_max_end_pos_list[rname]) - 1) 
    while((index >= 0) and 
          (gene_max_end_pos_list[rname][index] >= tstart)):

        if (min(tend, gene_end_pos_list[rname][index]) - 
            max(tstart, gene_start_pos_list[rname][index]) > 0):
            gnames.add(gname_list[rname][index])
        index -= 1

    return gnames
        
    
    
### Main
#############
[gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list] = parse_ref_file(ref_filename)

#chr = 'chr17'
#index =  gname_list[chr].index('USP32')
#print index
#print gene_end_pos_list[chr][index]
#print gene_start_pos_list[chr][index]
#print gene_max_end_pos_list[chr][index]
gpd_file = open(gpd_filename,'r')
output_file = open(output_filename ,'w')
temp_list = []



line_1 = gpd_file.readline().strip()
if (line_1):
    readname_1 = re.split(r'\+|\-', line_1.split()[0].split("|")[1])[0]

while (True):
    line_2 = gpd_file.readline().strip()
    if not line_2:
        break
    readname_2 = re.split(r'\+|\-', line_2.split()[0].split("|")[1])[0]
    if (readname_1 != readname_2):
        line_1 = line_2
        readname_1  = readname_2
        continue

    for line in [line_1, line_2]:
        line_ls = line.split('\t')
        gnames = get_gnames(line_ls[2], int(line_ls[4]), int(line_ls[5]),
                            gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list)
        if (len(gnames) == 0):
            gnames.add("NONE")
        
        output_file.write(line_ls[0])
        output_file.write('\t' + ','.join(list(gnames)) + '\t' + line_ls[2] + '\t')
        
    output_file.write('\n')
    # Read next line
    line_1 = gpd_file.readline().strip()
    if (line_1):
        readname_1 = re.split(r'\+|\-', line_1.split()[0].split("|")[1])[0]


gpd_file.close()
output_file.close()
################################################################################


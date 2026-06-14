#!/usr/bin/python

import sys
import os
import bisect
import psl_basics

if len(sys.argv) >= 5:
    psl_filename = sys.argv[1]
    skip_Nine = int(sys.argv[2])
    L_junction_limit = int(sys.argv[3])
    fusion_segment_len_threshold = int(sys.argv[4])
    fusion_overlap_threshold = int(sys.argv[5])
    fusion_gap_threshold = int(sys.argv[6])
    ref_filename = sys.argv[7]
    fusion_filename = sys.argv[8]
else:
    print("usage:blat_best.py psl_file skip_Nine L_junction_limit" +
          " fusion_segment_len_threshold fusion_overlap_threshold fusion_gap_threshold")
    sys.exit(1)

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
        
    
    
    
################################################################################
def process_temp_list(temp_list,
                      gname_list, gene_start_pos_list, 
                      gene_end_pos_list, gene_max_end_pos_list):
    ref_stat = 0
    
    # Select best ref_ls from single alignemt
    for result_ls in temp_list:
        stat = float(result_ls[0])/float(result_ls[10])
        if stat > ref_stat:
            ref_stat = stat
            ref_str = '\t'.join(result_ls)

    fusion_flag = False
    # Select best ref_ls from two alignemts
    # qstart/end index: 11/12
    # tstart/end index: 15/16
    # tname index: 13
    temp_list_len = len(temp_list)
    for idx in range(temp_list_len):
        # Check numbe rof matches
        gnames_1 = get_gnames(temp_list[idx][13], int(temp_list[idx][15]), int(temp_list[idx][16]),
                                 gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list)

        if (int(temp_list[idx][0]) < fusion_segment_len_threshold):
            continue
        line1 = "\t".join(temp_list[idx])
        entry1 = psl_basics.read_psl_entry(line1)
        for idx_2 in range(idx + 1, temp_list_len):
            # Check number of matches
            if (int(temp_list[idx_2][0]) < fusion_segment_len_threshold):
                continue
            if (temp_list[idx][13] == temp_list[idx_2][13]):
                #print ">" + str(temp_list)
                
                gnames_2 = get_gnames(temp_list[idx_2][13], int(temp_list[idx_2][15]), int(temp_list[idx_2][16]),
                                 gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list)
                # Check if they share same gene locus
                if ((len(gnames_1) > 0) and (len(gnames_2) > 0)):
                    if (len(gnames_1 & gnames_2) != 0):
                        continue
                else:
                    tstart = max(int(temp_list[idx][15]), int(temp_list[idx_2][15]))
                    tend = min(int(temp_list[idx][16]), int(temp_list[idx_2][16]))
                    # Note: it is expected intron shorter than L_junction_limit are not splited by blat
                    if ((tstart - tend) < L_junction_limit):
                        continue
            # Note: query coordinates have special handling in neg strands
            # Check psl format on genome USCSC website
            line2 = "\t".join(temp_list[idx_2])
            entry2 = psl_basics.read_psl_entry(line2)
            baseoverlap = psl_basics.query_coordinates_base_overlap_size(entry1,entry2)
            basegap = psl_basics.query_coordinates_gap_size(entry1,entry2)
            if (baseoverlap):
                if (baseoverlap <= fusion_overlap_threshold):
                    qlen = sum(entry1['blockSizes']) - baseoverlap
                    qlen_2 = sum(entry2['blockSizes']) - baseoverlap
                    if ((qlen < fusion_segment_len_threshold) or
                        (qlen_2 < fusion_segment_len_threshold)):
                        continue
                else:
                    continue
            elif (basegap > fusion_gap_threshold):
                continue
            stat = (float(temp_list[idx][0]) + float(temp_list[idx_2][0]))/float(temp_list[idx][10])
            if stat > ref_stat:
                fusion_psl.write('\t'.join(temp_list[idx]) + '\n')
                fusion_psl.write('\t'.join(temp_list[idx_2]) + '\n')
                fusion_flag = True


    if (fusion_flag):
        fusion_psl_single.write(ref_str + '\n')
        return ""
    else:
     return (ref_str)

################################################################################
[gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list] = parse_ref_file(ref_filename)

#chr = 'chr17'
#index =  gname_list[chr].index('USP32')
#print index
#print gene_end_pos_list[chr][index]
#print gene_start_pos_list[chr][index]
#print gene_max_end_pos_list[chr][index]


psl = open(psl_filename,'r')
fusion_psl = open(fusion_filename + '_pair','w')
fusion_psl_single = open(fusion_filename + '_single','w')
temp_list = []
Qname=""
i=0
for line in psl:
    if i < skip_Nine:
        i+=1
        continue
    ls = line.strip().split('\t')
    if (len(ls) != 21):
        sys.stderr.write("Warning: Invalid psl line in " + psl_filename + ": " + line)
        continue
    if Qname == ls[9]:
        temp_list.append(ls)
    else:
        if not Qname =="":
            result_str = process_temp_list(temp_list, gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list)
            if (result_str != ""):
                print result_str
        Qname = ls[9]
        temp_list = [ls]
    i+=1
result_str = process_temp_list(temp_list, gname_list, gene_start_pos_list, gene_end_pos_list, gene_max_end_pos_list)
if (result_str != ""):
    print result_str

fusion_psl.close()
fusion_psl_single.close()
################################################################################


#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2 :
    psl_filename = sys.argv[1]

else:
    print("usage: psl_filename")
    print("or ")
    sys.exit(1)
################################################################################

def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    return path, filename

psl = open(psl_filename,'r')
#junction_rd = open(junction_filename,'r')
#junction_wr = open(temp_foldername + "newname4_" + _junction_filename,'w')
ref_name_ls = []
ref_num_ls = []
read_name_ls = []
id = 1
left_seg = True
for line in psl:
    ls = line.strip().split("\t")
    identity = str(round(float(ls[0])/float(ls[10]),4))
    name = identity + "_" + ls[10]
    if (left_seg):
        #junc_line = junction_rd.readline().strip().split()
        #junction_wr.write('F' + str(id) +'\t' + '\t'.join(junc_line[1:]) + '\n')
        #ref_num_ls.append(int((junc_line[0].split('ref_'))[1]))
        #ref_name_ls.append(junc_line[0])
        #if (junc_line[2] == '+'):
        if (ls[8] == '+'):
            fusion_dir = '-'
            fusion_pos = ls[16]   #genome end point
        else:
            fusion_dir = '+'
            fusion_pos = ls[15]   #genome start point
        name += '|F' + str(id) + fusion_dir + fusion_pos + '|'
        left_seg = False
    else:
        #if (junc_line[2] == '-'):
        if (ls[8] == '+'):
            fusion_dir = '+'
            fusion_pos = ls[15]   #genome start point
        else:
            fusion_dir = '-'
            fusion_pos = ls[16]   #genome end point
        name += '|F' + str(id) + fusion_dir + fusion_pos + '|'
        left_seg = True
        id += 1
    if ls[9][-3:] == "ccs":
        name = name + "ccs" 
    ls[9] = name
    read_name_ls.append(name)
    print "\t".join(ls) 
psl.close()
#junction_rd.close()
#junction_wr.close()

"""
# Update pseudo-ref information
temp_foldername, _pseudo_ref_info_filename = GetPathAndName(pseudo_ref_info_filename) 
ref_num_ls.append(-1)  # Add invalid number to the end
file_rd = open(pseudo_ref_info_filename,'r')
file_wr = open(temp_foldername + "newname4_" + _pseudo_ref_info_filename,'w')
ref_id = 1
line_num = 1
for line in file_rd:
    line_ls = line.strip().split()
    if (ref_num_ls[0] == line_num):
        file_wr.write('\t'.join(['F' + str(ref_id), ref_name_ls[0], read_name_ls[0], read_name_ls[1]]) + '\t' + '\t'.join(line_ls) + '\n')
        del ref_num_ls[0]
        del ref_name_ls[0]
        del read_name_ls[0:2]
        ref_id += 1
    line_num +=1
file_rd.close()
file_wr.close()
"""
################################################################################
    

#!/usr/bin/python

import sys
import os
import re


if len(sys.argv) == 2:
    fusion_gdp_filename = sys.argv[1]
    fusion_single_filename = sys.argv[2]
    
else:
    sys.exit(1)
################################################################################
fusion_gdp_pair_file = open(fusion_gdp_filename, 'r')
#Note F# is in order with single file, so use that to keep single alignment lines that are filtered out
valid_fusion_lines = []
for line in fusion_gdp_pair_file:
    readname =  int(re.split(r'\+|\-', (line_1.split()[0].split('|')[1]))[0][1:])
    valid_fusion_lines.append(readname)
fusion_gdp_pair_file.close()

fusion_gdp_single_file_in = open(fusion_single_filename , 'r')
fusion_gdp_single_file_out = open(fusion_gdp_filename + '_single' , 'w')
line_num = 1
for line in fusion_gdp_single_file_in:
    if (line_num not in valid_fusion_lines):
        fusion_gdp_pair_file_out.write(line)
fusion_gdp_single_file_in.close()
fusion_gdp_single_file_out.close()

        
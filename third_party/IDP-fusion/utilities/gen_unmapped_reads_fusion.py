#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2:
    output_filename = sys.argv[1]
else:
    sys.exit(1)
################################################################################
""" Note: Expected SR read-name sorted sam file"""
mapped_reads = set()
for filename in sys.argv[2:-1]:
    print filename
    line_num = 0
    file = open(filename, 'r')
    for line in file:
        line_num += 1
        if (line_num % 1000000 == 0):
            print "line-num: " + str(line_num)
        if (line[0] == '@'):
            continue
        
        fields = line.split('\t')
        read_name = fields[0]
        if (fields[5] != '*'):
            mapped_reads.add(read_name)
    file.close()
if (fields[5] != '*'):
    mapped_reads.add(read_name)

output_file = open(output_filename, 'w')
# Parse the last sam file and store the unmapped reads
filename = sys.argv[-1]
print filename
line_num = 0
file = open(filename, 'r')
for line in file:
    line_num += 1
    if (line_num % 1000000 == 0):
        print "line-num: " + str(line_num)
    if (line[0] == '@'):
        continue
        
    fields = line.split('\t')
    read_name = fields[0]
    read_seq = fields[9]

    if ((fields[5] == '*') and
        (read_name not in mapped_reads)):
        output_file.write('>' + read_name + '\n')
        output_file.write(read_seq + '\n')

if (fields[5] == '*' and
    (read_name not in mapped_reads)):
    output_file.write('>' + read_name + '\n')
    output_file.write(read_seq + '\n')
file.close()
output_file.close()

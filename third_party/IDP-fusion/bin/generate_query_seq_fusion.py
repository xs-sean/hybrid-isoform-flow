#!/usr/bin/python
import sys
import os


if len(sys.argv) >= 4:
    read_filename = sys.argv[1]
    fusion_psl_filename = sys.argv[2]
    output_filename = sys.argv[3]
else:
    sys.exit(1)


psl_file = open(fusion_psl_filename, 'r')
# Get list of read-names
LR_names_dict = dict()
for line in psl_file:
    LR_names_dict[line.split()[9]] = ''
psl_file.close()

########
file = open(read_filename, 'r')
for line in file:
    line = line.strip()
    if (line[0] == '>'):
        if (LR_names_dict.has_key(line[1:]) and
            (LR_names_dict[line[1:]] == '')):
            LR_name = line[1:]
        else:
            # There could be duplicate names (use sequence from first occurrence)
            LR_name = ''
    else:
        if (LR_name != ''):
            LR_names_dict[LR_name] = LR_names_dict[LR_name] + line
file.close()

########
psl_file = open(fusion_psl_filename, 'r')
output_file = open(output_filename, 'w')
printed_reads_ls = []
for line in psl_file:

    fields = line.split('\t')
    # Add margin to ends
    margin = 0
    read_name = fields[9]
    seq = LR_names_dict[read_name][max(0, int(fields[11]) - margin):min(int(fields[10]), int(fields[12]) + margin)]
    read_name = read_name + '/' + fields[11] + '_' + fields[12]
    if (read_name not in printed_reads_ls):
        output_file.write('>' + read_name  + '\n')
        output_file.write(seq + '\n')
        printed_reads_ls.append(read_name)

output_file.close()
psl_file.close()

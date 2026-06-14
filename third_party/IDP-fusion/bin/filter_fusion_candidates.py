#!/usr/bin/python
import sys
import os


if len(sys.argv) >= 4:
    psl_filename = sys.argv[1]
    fusion_psl_filename = sys.argv[2]
    min_similarity_gap_perc =  float(sys.argv[3]) / 100
else:
    sys.exit(1)

# Consider some margin to check best alignment is the original one in fusion psl
margine = 2

################
psl_file = open(psl_filename, 'r')
# Get list of valid segments from re-mapping (value: best-alignment line)
readnames_dict = dict()
stat_ls = []
best_stat = 0
prev_read_name = ''
best_alignemnt = ''
read_name = ''
for line in psl_file:
    fields = line.split()
    read_name = fields[9] 
    if (prev_read_name != read_name):
        if (prev_read_name != ''):
            stat_ls = sorted(stat_ls)
            if (len(stat_ls) > 1):
                # Check if it is the best unique alignment
                if ((stat_ls[-1] - stat_ls[-2]) > min_similarity_gap_perc):
                    readnames_dict[prev_read_name] = best_alignemnt
                else:
                    readnames_dict[prev_read_name] = False
            else:
                readnames_dict[prev_read_name] = best_alignemnt
                
        stat_ls = []
        prev_read_name = read_name
        best_stat = 0
        best_alignemnt = ''
    stat = float(fields[0])/float(fields[10])
    stat_ls.append(stat)
    if (stat > best_stat):
        best_alignemnt = line
        best_stat = stat
stat_ls = sorted(stat_ls)
if (len(stat_ls) > 1):
    if ((stat_ls[-1] - stat_ls[-2]) > min_similarity_gap_perc):
        # Save the best alignemt
        readnames_dict[read_name] = best_alignemnt
    else:
        readnames_dict[read_name] = False
else:
    readnames_dict[read_name] = best_alignemnt
psl_file.close()

# Filter out alignment lines in pair psl file
fusion_readnames = []
psl_file_pair = open(fusion_psl_filename + '_pair', 'r')
output_psl_file_pair = open(fusion_psl_filename + '_pair_filtered', 'w')
prev_read_name = ''
best_stat = 0
while (True):
    line_f1 = psl_file_pair.readline().strip()
    if not line_f1:
        break
    line_f2 = psl_file_pair.readline().strip()
    fields = [line_f1.split(), line_f2.split()]
    # Check read-names consistency
    if (fields[0][9] != fields[1][9]):
        print "Err: Invalid fusion_*.psl file format"
        exit(1)

    if (prev_read_name != fields[0][9]):
        if (prev_read_name != ''):
            if (best_stat > 0):
                output_psl_file_pair.write(line_1 + '\n')
                output_psl_file_pair.write(line_2 + '\n')
                fusion_readnames.append(prev_read_name)
                
        prev_read_name = fields[0][9]
        best_stat = -1
        line_1 = ''
        line_2 = ''
        
    valid_flag = True
    # Check if the paired alignments are valid fusion-candidates
    for i in [0, 1]:
        readname = fields[i][9] + '/' + fields[i][11] + '_' + fields[i][12]
        if ((not readnames_dict.has_key(readname)) or (readnames_dict[readname] == False)):
            valid_flag = False
        else:
            
            alignment_fields = readnames_dict[readname].split()
            if (abs(int(alignment_fields[0]) - int(fields[i][0])) > margine):
                valid_flag = False

            # chr-name 
            if (alignment_fields[13] != fields[i][13]):
                valid_flag = False
            # chr start-pos
            if (abs(int(alignment_fields[15]) - int(fields[i][15])) > margine):
                valid_flag = False
    if (valid_flag):
        stat = (float(fields[0][0]) + float(fields[1][0])) / float(fields[0][10])
    else:
        stat = -1
    if (stat > best_stat):
        # Save them in read start-pos order
        if (int(fields[0][11]) < int(fields[1][11])):
            line_1 = line_f1
            line_2 = line_f2
        else:
            line_1 = line_f2
            line_2 = line_f1
        best_stat = stat
if (best_stat > 0):
    output_psl_file_pair.write(line_1 + '\n')
    output_psl_file_pair.write(line_2 + '\n')

psl_file_pair.close()
output_psl_file_pair.close()


# Filter out alignment lines in pair psl file
psl_file_single = open(fusion_psl_filename + '_single', 'r')
output_psl_file_single = open(fusion_psl_filename + '_single_filtered', 'w')
for line in psl_file_single:
    fields = line.split('\t')
    if (fields[9] not in fusion_readnames):
        output_psl_file_single.write(line)
psl_file_single.close()
output_psl_file_single.close()
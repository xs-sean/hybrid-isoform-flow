#!/usr/bin/python
import sys
import re
import string
#from constants import *
import re


valid_cigar = set("0123456789MNID")
read_len_margin = 0

### 
##########
def compute_min_ljust(str_value):
    return max(20, (len(str_value)/10 + 1) * 10)


### Due to redundancy, isoforms might not be generated for a fusion segment, 
### Need to find the super/master LR 
###########
def find_super_fusion_id(fusion_name, read_mapping_dict, valid_readname_ls):
    key_2 = fusion_name
    key_1 = "invalid"
    while (key_1 != key_2):
        key_1 = key_2
        if (read_mapping_dict.has_key(key_1)):
            key_2 = read_mapping_dict[key_1]
        else:
            break
    if (key_1 in valid_readname_ls):
        return key_1
    else:
        return -1

###
###########
def parse_fusion_genenames(fusion_compatible_filename, fusion_gpd_filename, readmapping_filename):
    
    fusion_isofroms_pair_ls = []
    MLE_input_dict = {}

    file = open(readmapping_filename, 'r')
    read_mapping_dict = dict()
    for line in file:
        fields = line.split()
        read_mapping_dict[fields[0]] = fields[1]
    file.close()
    

    file = open(fusion_gpd_filename, 'r')
    read_mate_dict = dict()
    readnames_ls = []
    while True:
        line_1_ls = file.readline().strip().split()
        line_2_ls = file.readline().strip().split()
        if line_1_ls == []:
            break
        read_mate_dict[line_2_ls[0]] = line_1_ls[0]
        read_mate_dict[line_1_ls[0]] = line_2_ls[0]
        readnames_ls += [line_1_ls[0], line_2_ls[0]]
    file.close()


    fusion_gpd_dict = dict()
    fusion_gpd_file = open(fusion_compatible_filename, 'r')
    for line in fusion_gpd_file:
        fields = line.split()
        fusion_gpd_dict[fields[0]] = line
        fields = line.split()
        if (not read_mapping_dict.has_key(fields[0])):
            read_mapping_dict[fields[0]] = fields[0]
    fusion_gpd_file.close()
    
    processed_reads_ls = []  # To avoid duplicates
    valid_readname_ls = fusion_gpd_dict.keys()
    compreadname2readname_mapping = {}
    for readname_origin in readnames_ls:
        readname = find_super_fusion_id(readname_origin, read_mapping_dict, valid_readname_ls)
        if readname == -1:
            continue
        fields = fusion_gpd_dict[readname].split()
        [fusion_id, fusion_dir, fusion_pos]= re.split(r"(\+|\-)", (fields[0].split("|"))[1])
        if (fusion_dir == '+'):   # Change it to 1-based
            fusion_pos = str(1 + int(fusion_pos))
        MLE_input_dict[fields[1]] = []
        
        paired_readname = find_super_fusion_id(read_mate_dict[readname_origin], read_mapping_dict, valid_readname_ls)

        # Keep mapping of original readnames to uniqueue ones (for reference/debugging)
        if (paired_readname != -1):
            if (not compreadname2readname_mapping.has_key(readname)):
                compreadname2readname_mapping[readname] = {}
            if (not compreadname2readname_mapping[readname].has_key(paired_readname)):
                compreadname2readname_mapping[readname][paired_readname] = set()
            compreadname2readname_mapping[readname][paired_readname].add(re.split(r"(\+|\-)", (readname_origin.split("|"))[1])[0])

            if (not compreadname2readname_mapping.has_key(paired_readname)):
                compreadname2readname_mapping[paired_readname] = {}
            if (not compreadname2readname_mapping[paired_readname].has_key(readname)):
                compreadname2readname_mapping[paired_readname][readname] = set()
            compreadname2readname_mapping[paired_readname][readname].add(re.split(r"(\+|\-)", (readname_origin.split("|"))[1])[0])


        if (((str(paired_readname)+'_'+readname) not in  processed_reads_ls) and
            ((readname+'_'+str(paired_readname)) not in  processed_reads_ls)):
            
            # gene-names and isoform-points will be added
            append_ls = [[fields[0], fusion_dir, fusion_pos, fields[1], 0, 0]]
            if (fusion_dir == '+'):
                points_temp =  fields[9].split(',')[:-1]
                append_ls[0][4] = [(int(i)+1) for i in points_temp]
                points_temp =  fields[10].split(',')[:-1]
                append_ls[0][5] = [int(i) for i in points_temp]
            elif (fusion_dir == '-'):
                points_temp =  fields[9].split(',')[:-1]
                append_ls[0][5] = list(reversed([(int(i)+1) for i in points_temp]))
                points_temp =  fields[10].split(',')[:-1]
                append_ls[0][4] = list(reversed([int(i) for i in points_temp]))
            append_ls[0][4][0] = int(fusion_pos)
        
            if (paired_readname != -1):
                
                processed_reads_ls.append(readname+'_'+str(paired_readname))
                paired_fields = fusion_gpd_dict[paired_readname].split()
                
                MLE_input_dict[paired_fields[1]] = []
                [paired_fusion_id, paired_fusion_dir, paired_fusion_pos]= re.split(r"(\+|\-)", (paired_fields[0].split("|"))[1])
                if (paired_fusion_dir == '+'):   # Change it to 1-based
                    paired_fusion_pos = str(1 + int(paired_fusion_pos))
                append_ls.append([paired_fields[0], paired_fusion_dir, paired_fusion_pos, paired_fields[1], 0, 0])
                if (paired_fusion_dir == '+'):
                    points_temp =  paired_fields[9].split(',')[:-1]
                    append_ls[1][4] = [(int(i)+1) for i in points_temp]
                    points_temp =  paired_fields[10].split(',')[:-1]
                    append_ls[1][5] = [int(i) for i in points_temp]
                elif (paired_fusion_dir == '-'):
                    points_temp =  paired_fields[9].split(',')[:-1]
                    append_ls[1][5] = list(reversed([(int(i)+1) for i in points_temp]))
                    points_temp =  paired_fields[10].split(',')[:-1]
                    append_ls[1][4] = list(reversed([int(i) for i in points_temp]))
                append_ls[1][4][0] = int(paired_fusion_pos)
                
            else:
                pass
                #append_ls.append([0, 0, 0, 0, 0, 0])
        
            fusion_isofroms_pair_ls.append(append_ls)
            
    return [fusion_isofroms_pair_ls, MLE_input_dict, compreadname2readname_mapping]
    
    
### Store fusion gene information from MLE input file
#########
def get_fusion_genomes(input_filename, MLE_input_dict):
    
        
    input_file = open(input_filename, 'r')
    line_index = 0;
    input_file.readline()   # Pass the top line
    wr_flag = False
    for line in input_file:
        fields = line.split()
        if (line_index == 0):
            gname = fields[0]
            num_isoforms = int(fields[1])
            rname = fields[2]
            if (gname in MLE_input_dict.keys()):
                wr_flag = True
        line_index += 1
        if (wr_flag):
            MLE_input_dict[gname].append(line.strip())
        if (line_index == (9 + num_isoforms)):
            line_index = 0
            wr_flag = False
            
    input_file.close()
    
    for gname in MLE_input_dict.keys():
        if MLE_input_dict[gname] == []:
            del MLE_input_dict[gname]
    return MLE_input_dict


"""
###
### Note: fusion_isofroms_pair_ls elements that have '0' gnames, do not have generated isoforms!
###########
def map_id2genome(fusion_compatible_filename,
                  fusion_isofroms_pair_ls, MLE_input_dict):
    
    # Remove gene-names that are not generated in isoform_construction step
    for gname in MLE_input_dict.keys():
        if (len( MLE_input_dict[gname]) < 5):
            del MLE_input_dict[gname]
            continue
    
    file = open(fusion_compatible_filename, 'r')
    for line in file:
        fields = line.split()
        [fusion_id, fusion_dir, fusion_pos]= re.split(r"(\+|\-)", (fields[0].split("|"))[1])
        if (fusion_isofroms_pair_ls.has_key(fusion_id)):
            if (fusion_isofroms_pair_ls[fusion_id][0][0] == fields[0]):
                read_idx = 0
            elif (fusion_isofroms_pair_ls[fusion_id][1][0] == fields[0]):
                read_idx = 1
                            
            fusion_isofroms_pair_ls[fusion_id][read_idx][3] = fields[1]   #gene-name
            
            if (fusion_isofroms_pair_ls[fusion_id][read_idx][1] == '+'):
                points_temp =  fields[9].split(',')[:-1]
                fusion_isofroms_pair_ls[fusion_id][read_idx][4] = [(int(i)+1) for i in points_temp]
                points_temp =  fields[10].split(',')[:-1]
                fusion_isofroms_pair_ls[fusion_id][read_idx][5] = [int(i) for i in points_temp]
            elif (fusion_isofroms_pair_ls[fusion_id][read_idx][1] == '-'):
                points_temp =  fields[9].split(',')[:-1]
                fusion_isofroms_pair_ls[fusion_id][read_idx][5] = list(reversed([(int(i)+1) for i in points_temp]))
                points_temp =  fields[10].split(',')[:-1]
                fusion_isofroms_pair_ls[fusion_id][read_idx][4] = list(reversed([int(i) for i in points_temp]))
            fusion_isofroms_pair_ls[fusion_id][read_idx][4][0] = int(fusion_isofroms_pair_ls[fusion_id][read_idx][2])

        read_mapping_dict[fields[0]] = fields[0]
    file.close()
    

    file = open(readmapping_filename, 'r')
    for line in file:
        fields = line.split()
        read_mapping_dict[fields[0]] = fields[1]
    file.close()
    
    for fusion_id in fusion_isofroms_pair_ls.keys():
        
        for read_idx in [0, 1]:
            if fusion_isofroms_pair_ls[fusion_id][read_idx][3] == 0:    # gene-name is not assigned
                super_fusion = find_super_fusion_id(fusion_isofroms_pair_ls[fusion_id][read_idx][0], read_mapping_dict)
                if (super_fusion == -1):
                    break  # Could not find associated gname
                super_fusion_id = re.split(r"(\+|\-)", (super_fusion.split("|"))[1])[0]
                if (fusion_isofroms_pair_ls[super_fusion_id][0][0] == super_fusion):
                    super_read_idx = 0
                elif (fusion_isofroms_pair_ls[super_fusion_id][1][0] == super_fusion):
                    super_read_idx = 1
                fusion_isofroms_pair_ls[fusion_id][read_idx][3] = fusion_isofroms_pair_ls[super_fusion_id][super_read_idx][3]
                fusion_isofroms_pair_ls[fusion_id][read_idx][4] = fusion_isofroms_pair_ls[super_fusion_id][super_read_idx][4]
                fusion_isofroms_pair_ls[fusion_id][read_idx][5] = fusion_isofroms_pair_ls[super_fusion_id][super_read_idx][5]
                
    valid_genenames = MLE_input_dict.keys()
    for fusion_id in fusion_isofroms_pair_ls.keys():
        valid_flag = True
        if (fusion_isofroms_pair_ls[fusion_id][0][3] not in valid_genenames):
            valid_flag = False
        elif (fusion_isofroms_pair_ls[fusion_id][1][3] not in valid_genenames):
            valid_flag = False
              
        if (not valid_flag):
            del fusion_isofroms_pair_ls[fusion_id]
                
    
    
    return fusion_isofroms_pair_ls
"""

### Merge genes data which are part of a fusion gene 
###########
def merge_fusion_genes(MLE_input_dict, fusion_isofroms_pair_ls):
    fusion_MLE_input_dict = dict()
    fusion_gnames_map = dict()
    
    for gname in MLE_input_dict.keys():
        fusion_gnames_map[gname] = gname
        
    # Merge genomes part of fusion genes
    for fusion_isofroms_pair in fusion_isofroms_pair_ls:
        
        if (not fusion_gnames_map.has_key(fusion_isofroms_pair[0][3])):
                gname_1 = ''
        else:
            gname_1 = fusion_gnames_map[fusion_isofroms_pair[0][3]]
        if len(fusion_isofroms_pair) == 1:
            # not paired
            gname_2 = ''
        else:
            if (not fusion_gnames_map.has_key(fusion_isofroms_pair[1][3])):
                gname_2 = ''
            else:
                gname_2 = fusion_gnames_map[fusion_isofroms_pair[1][3]]

        if (gname_1 == gname_2):
            continue # They are already merged 
        
        # Check if gname is generated as part of isoform construction (e.g. not skipped because num isoforms > limit)
        fusion_gname = ''
        if ((not fusion_MLE_input_dict.has_key(gname_1)) and
            (not MLE_input_dict.has_key(gname_1))):
            gname_1 = ''
        else:
            fusion_gname += gname_1 + '/'
            
        if ((not fusion_MLE_input_dict.has_key(gname_2)) and
            (not MLE_input_dict.has_key(gname_2))):
            gname_2 = ''
        else:
            fusion_gname += gname_2 + '/'
        
        if (fusion_gname == ''):
            continue
        else:
            fusion_gname = fusion_gname[:-1]
        
       
        gname_ls = list(set(fusion_gname.split('/')))  # Remove redundancy
        fusion_gname = '/'.join(gname_ls)
        for gname in gname_ls:
            fusion_gnames_map[gname] = fusion_gname
        
        # Reconstructs the joint genome
        if (fusion_MLE_input_dict.has_key(gname_1)):  
            del fusion_MLE_input_dict[gname_1]
        if (fusion_MLE_input_dict.has_key(gname_2)):
            del fusion_MLE_input_dict[gname_2]

        num_isoforms = 0
        num_isoforms_ls = []
        chr_name = ''
        for gname in gname_ls:
            num_isoforms += int(MLE_input_dict[gname][0].split()[1])
            num_isoforms_ls.append(int(MLE_input_dict[gname][0].split()[1]))
            chr_name += MLE_input_dict[gname][0].split()[2] + ','
        chr_name = chr_name[:-1]  # Exclude last ','
        # gname, num_isoforms, chr
        line = '\t'.join([fusion_gname, ','.join([str(num_isoform) for num_isoform in num_isoforms_ls]), chr_name])
        fusion_MLE_input_dict[fusion_gname] = [line]
        
        for i in range(1, 4):
            line = ''
            for gname in gname_ls:
                line += MLE_input_dict[gname][i] + '\t'
            fusion_MLE_input_dict[fusion_gname].append(line.strip())
        
        line = ''
        i = 0
        for gname in gname_ls:
            line += MLE_input_dict[gname][4].replace('P', string.uppercase[i]) + '\t'
            i += 1
        fusion_MLE_input_dict[fusion_gname].append(line.strip())
        
        for i in [5]:
            line = ''
            for gname in gname_ls:
                line += MLE_input_dict[gname][i] + '\t'
            fusion_MLE_input_dict[fusion_gname].append(line.strip())
        
        line = ''
        i = 0
        for gname in gname_ls:
            line += MLE_input_dict[gname][6].replace('P', string.uppercase[i]) + '\t'
            i += 1
        fusion_MLE_input_dict[fusion_gname].append(line.strip())
        
        accum_num_regions = 0
        total_num_regions = len(line.split())
        for idx_1 in range(len(gname_ls)):
            num_regions = len(MLE_input_dict[gname_ls[idx_1]][6].split())
            for idx_2 in range(num_isoforms_ls[idx_1]):
                line = '\t'.join([str(0)] * accum_num_regions) + '\t'
                line += MLE_input_dict[gname_ls[idx_1]][7 + idx_2] + '\t'
                line += '\t'.join([str(0)] * (total_num_regions - accum_num_regions - num_regions))
                fusion_MLE_input_dict[fusion_gname].append(line.strip())
            accum_num_regions += num_regions
        
        for i in [-2, -1]:
            line = ''
            for gname in gname_ls:
                line += MLE_input_dict[gname][i] + '\t'
            fusion_MLE_input_dict[fusion_gname].append(line.strip())
        
    # Remove fusion_gnames_map that are not in fusion_MLE_input_dict
    for gname in fusion_gnames_map.keys():
        fusion_gname = fusion_gnames_map[gname]
        if (not fusion_MLE_input_dict.has_key(fusion_gname)):
            del fusion_gnames_map[gname]
    return [fusion_MLE_input_dict, fusion_gnames_map]

### Output valid regions around fusion points
### Input: start position index of  each within fusion gene,
###        fusion genes information
###        MLE-format information for fusion gene
##########
def get_fusion_regions(genes_start_point_idx, fusion_id_info, 
                       fusion_genename_info):
    
    # Keep points in order to avoid duplicate regions
    if (genes_start_point_idx[0][0] < genes_start_point_idx[1][0]):
        read_idx_ls = [0, 1]
    else:
        read_idx_ls = [1, 0]
    fusion_points_start = [0, 0]
    fusion_points_end = [0, 0]
    # Get exon boundaries of fusion genes
    for read_idx in read_idx_ls:
        fusion_points_start[read_idx] = fusion_id_info[read_idx][4]
        fusion_points_end[read_idx] = fusion_id_info[read_idx][5]
    
    fusion_regions_ls = []
    fusion_regions_length_ls = []
    point_ls = [int(i) for i in fusion_genename_info[5].split()]
    point_names_ls = fusion_genename_info[4].split()
    
    region_dict = {}
    segment_len = [(READ_LEN - READ_JUNC_MIN_MAP_LEN), READ_JUNC_MIN_MAP_LEN]
    while (True):
        
        subregion = [[], []]
        for read_idx in read_idx_ls:
            is_valid_region = True
            # First extend left gene by READ_LEN - READ_JUNC_MIN_MAP_LEN   
            idx = -1
            len_temp = 0
            while (len_temp < segment_len[read_idx]):
                 idx += 1
                 prev_len_temp = len_temp
                 len_temp += abs(fusion_points_start[read_idx][idx] - fusion_points_end[read_idx][idx]) + 1  # start and end points are inclusive
                 if ((idx + 1) == len(fusion_points_start[read_idx])):
                     break
            # Update seg-length (in case, initially it is less than (READ_LEN - READ_JUNC_MIN_MAP_LEN))
            segment_len[read_idx] = min(len_temp, segment_len[read_idx])
            # Last exon segment should have minimum READ_JUNC_MIN_MAP_LEN length to have reliable read mapping
            if ((segment_len[read_idx] - prev_len_temp) < READ_JUNC_MIN_MAP_LEN):
                is_valid_region = False
                break
            if (fusion_id_info[read_idx][1] == '+'):
                end_point = fusion_points_start[read_idx][idx] + (segment_len[read_idx] - prev_len_temp) - 1
            else:
                end_point = fusion_points_start[read_idx][idx] - (segment_len[read_idx] - prev_len_temp) + 1
            start_point = fusion_points_start[read_idx][0]
            idx = 0
            
            # Return empty dictionary if isorom for this fusio gene is not generated
            if (fusion_points_start[read_idx][0] not in point_ls[genes_start_point_idx[read_idx][0]: genes_start_point_idx[read_idx][1]]):
                return dict()
            
            point_idx = point_ls[genes_start_point_idx[read_idx][0]: genes_start_point_idx[read_idx][1]].index(start_point) + genes_start_point_idx[read_idx][0]
            subregion[read_idx].append(point_names_ls[point_idx])
            if (fusion_id_info[read_idx][1] == '+'):
                point_idx += 1
            else:
                point_idx -= 1
            start_point = point_ls[point_idx]
            if (fusion_id_info[read_idx][1] == '+'):
                while (start_point <= end_point):
                    if (point_ls[point_idx] >= fusion_points_start[read_idx][idx]):
                        if (point_ls[point_idx] <= fusion_points_end[read_idx][idx]):
                            subregion[read_idx].append(point_names_ls[point_idx])
                            point_idx += 1
                            start_point = point_ls[point_idx]
                        else:
                            idx += 1
                    else:
                        point_idx += 1
                        start_point = point_ls[point_idx]
            elif (fusion_id_info[read_idx][1] == '-'):
                while (start_point >= end_point):
                    if (point_ls[point_idx] >= fusion_points_end[read_idx][idx]):
                        if (point_ls[point_idx] <= fusion_points_start[read_idx][idx]):
                            subregion[read_idx].append(point_names_ls[point_idx])
                            point_idx -= 1
                            start_point = point_ls[point_idx]
                        else:
                            idx += 1
                    else:
                        point_idx -= 1
                        start_point = point_ls[point_idx]
                
        if (is_valid_region):
            region = '-'.join(subregion[read_idx_ls[0]][::-1]) + '-' + '-'.join(subregion[read_idx_ls[1]])
            if (region_dict.has_key(region)):
                region_dict[region] += 1
            else:
                region_dict[region] = 1
                
        segment_len[0] -= 1
        segment_len[1] += 1
        if (segment_len[0] < READ_JUNC_MIN_MAP_LEN):
            break

    return region_dict
    
    
    
###
###########
def add_fusion_regions(fusion_MLE_input_dict, fusion_gnames_map, 
                        fusion_isofroms_pair_ls):
    
    for fusion_isofroms_pair in fusion_isofroms_pair_ls:
        
        if (len(fusion_isofroms_pair) == 1):
            #Not paired
            continue
        
        try:
            gname_1 = fusion_isofroms_pair[0][3]
            gname_2 = fusion_isofroms_pair[1][3]
        
            fusion_gene_name_1 = fusion_gnames_map[gname_1]
            fusion_gene_name_2 = fusion_gnames_map[gname_2]
        except KeyError:
            continue
        
        if ( fusion_gene_name_1 != fusion_gene_name_2):
            print "Err: expected to have same fusion gene-name"
            exit(1)  
        
        fusion_gene_name = fusion_gene_name_1
        gene_names_ls = (fusion_MLE_input_dict[fusion_gene_name][0].split())[0].split('/')
        point_names_ls = fusion_MLE_input_dict[fusion_gene_name][4].split()
        genes_start_point_names = [0]  # point_name starting index for each gene
        alpha = point_names_ls[0][0]
        for idx in range(len(point_names_ls)):
            if (point_names_ls[idx][0] != alpha):  # Note: we might not have P0 if exon is of length 1 (only odd point is present then)
                genes_start_point_names.append(idx)
                alpha = point_names_ls[idx][0]

        if (len(genes_start_point_names) != len(gene_names_ls)):
            print "Err: Unmatched list lengths"
            exit(1)
        genes_start_point_names.append(len(point_names_ls))
        # Add this info to the fusion gene dict to be used later (if not already)
        if (len(fusion_MLE_input_dict[fusion_gene_name][0].split()) < 4):
            fusion_MLE_input_dict[fusion_gene_name][0] += '\t' + ','.join([str(i) for i in genes_start_point_names])
        
        genes_start_point_idx = [gene_names_ls.index(gname_1), gene_names_ls.index(gname_2)]
        genes_start_point_names = [[genes_start_point_names[genes_start_point_idx[0]], genes_start_point_names[genes_start_point_idx[0] +1]],
                                   [genes_start_point_names[genes_start_point_idx[1]], genes_start_point_names[genes_start_point_idx[1] +1]]]
        # Get list of fusion regions
        region_dict = get_fusion_regions(genes_start_point_names, fusion_isofroms_pair, 
                                         fusion_MLE_input_dict[fusion_gene_name])
        
        num_isoforms = [int(i) for i in ((fusion_MLE_input_dict[fusion_gene_name][0].split())[1]).split(',')]
        num_isoforms = sum(num_isoforms)
        regions_ls = fusion_MLE_input_dict[fusion_gene_name][6].split()
        regions_len_ls = fusion_MLE_input_dict[fusion_gene_name][7 + num_isoforms].split()
        fusion_regions = region_dict.keys()
        for fusion_region in fusion_regions:
            if (fusion_region in regions_ls):
                del region_dict[fusion_region]
            else:
                regions_ls.append(fusion_region)
                regions_len_ls.append(str(region_dict[fusion_region]))
        fusion_MLE_input_dict[fusion_gene_name][6] = '\t'.join(regions_ls)
        fusion_MLE_input_dict[fusion_gene_name][7 + num_isoforms] = '\t'.join(regions_len_ls)
        
        # Add this info to fusion pairs
        fusion_isofroms_pair.append(fusion_regions)
        
    return fusion_MLE_input_dict

 
###
### Removes isoforms used as part of fusion isoforms 
### Fill up isoform matrix for missing entries in fusion region for normal isoforms
##########
def finalize_fusion_genes(fusion_MLE_input_dict):
    
    fusion_gname_ls = fusion_MLE_input_dict.keys()
    for fusion_gname in fusion_gname_ls:
        line_ls = fusion_MLE_input_dict[fusion_gname][0].split()
        num_isforms = sum([int(i) for i in line_ls[1].split(',')])
        isoform_names_ls = fusion_MLE_input_dict[fusion_gname][1].strip().split()
        num_regions = len(fusion_MLE_input_dict[fusion_gname][6].split())

        # Remove isoforms ending in .f
        iso_idx = 0
        isoform_idx_to_remove = []
        for isoform_name in isoform_names_ls:
            if (isoform_name[-2:] == ".f"):
                isoform_idx_to_remove.append(iso_idx)
            iso_idx += 1
        isoform_idx_to_remove = list(reversed(isoform_idx_to_remove))
        isoform_names_ls = fusion_MLE_input_dict[fusion_gname][1].split()
        isoform_mark_ls = fusion_MLE_input_dict[fusion_gname][2].split()
        isoform_len_ls = fusion_MLE_input_dict[fusion_gname][3].split()
        for isoform_idx in isoform_idx_to_remove:
            del fusion_MLE_input_dict[fusion_gname][7 + isoform_idx]
            del isoform_names_ls[isoform_idx]
            del isoform_len_ls[isoform_idx]
            del isoform_mark_ls[isoform_idx]
            num_isforms -= 1
        fusion_MLE_input_dict[fusion_gname][1] = '\t'.join(isoform_names_ls)
        fusion_MLE_input_dict[fusion_gname][2] = '\t'.join(isoform_len_ls)
        # Mark fusion isoforms as known
        fusion_MLE_input_dict[fusion_gname][3] = '\t'.join(isoform_mark_ls + ['1'] * (num_isforms - len(isoform_mark_ls)))
        
        for isoform_idx in range(num_isforms):    # Add 0 for normal isforms in fusion-region indices
            isoform_ls = fusion_MLE_input_dict[fusion_gname][7 + isoform_idx].split()
            if (len(isoform_ls) < num_regions):
                isoform_ls += ['0'] * (num_regions - len(isoform_ls))
                fusion_MLE_input_dict[fusion_gname][7 + isoform_idx] = '\t'.join(isoform_ls)
            

        # Add 0 readcount for fusion regions
        read_count_ls = fusion_MLE_input_dict[fusion_gname][-1].strip().split()
        fusion_MLE_input_dict[fusion_gname][-1] = '\t'.join(read_count_ls + ['0'] * (num_regions - len(read_count_ls)))
        line_ls[1] = str(num_isforms)
        fusion_MLE_input_dict[fusion_gname][0] = '\t'.join(line_ls + [str(len(read_count_ls))])  #Inlude the start point for fusion regions
    
    return fusion_MLE_input_dict
        
    
###
###########
def add_fusion_isoforms(fusion_MLE_input_dict, fusion_gnames_map, 
                        fusion_isofroms_pair_ls, isoforms_readnames_filename):

    isoforms_readnames_file = open(isoforms_readnames_filename, 'r')
    fusion_isoforms_map = dict()
    for line in isoforms_readnames_file:
        fields = line.strip().split()
        if (fusion_isoforms_map.has_key(fields[1])):
            fusion_isoforms_map[fields[1]].append(fields[0])
        else:
            fusion_isoforms_map[fields[1]] = [fields[0]]
    isoforms_readnames_file.close()
    
    processed_fusion_isoforms = []
    for fusion_isofroms_pair in fusion_isofroms_pair_ls:

        gname = []
        gname_1 = fusion_isofroms_pair[0][3]
        if (fusion_gnames_map.has_key(gname_1)):
            fusion_gene_name_1 = fusion_gnames_map[gname_1]
            gname.append(gname_1) 
        else:
            fusion_gene_name_1 = ''
        fusion_gene_name_2 = ''
        if (len(fusion_isofroms_pair) > 1):
            gname_2 = fusion_isofroms_pair[1][3]
            if (fusion_gnames_map.has_key(gname_2)):
                fusion_gene_name_2 = fusion_gnames_map[gname_2]
                gname.append(gname_2)
        

        if ((fusion_gene_name_1 == fusion_gene_name_2) and
            (fusion_gene_name_1 != '')):
            fusion_gene_name = fusion_gene_name_1
            read_idx_ls = [0, 1]
        else:
            continue
            
            
        gene_names_ls = (fusion_MLE_input_dict[fusion_gene_name][0].split())[0].split('/')
        isoforms_names_ls = fusion_MLE_input_dict[fusion_gene_name][1].split()
        point_names_ls = fusion_MLE_input_dict[fusion_gene_name][4].split()
        points_ls = [int(i) for i in fusion_MLE_input_dict[fusion_gene_name][5].split()]
        gene_regions_ls = fusion_MLE_input_dict[fusion_gene_name][6].split()
        num_gene_points = len(points_ls)
        num_gene_regions = len(gene_regions_ls)
        gene_point_name_mapping = dict()
        for idx in range(num_gene_points):
            gene_point_name_mapping[point_names_ls[idx]] = points_ls[idx]
        
        num_isoforms_ls = ((fusion_MLE_input_dict[fusion_gene_name][0].split())[1]).split(',')
        num_isoforms_ls = [int(i) for i in num_isoforms_ls]
        num_isoforms = sum(num_isoforms_ls)
        
        split_ls = fusion_MLE_input_dict[fusion_gene_name][0].split()
        if (len(split_ls) >= 4):
            start_points_indices = [int(i) for i in ((split_ls)[3]).split(',')]
        else:
            # This not realy merged gene
            start_points_indices = [0, num_gene_regions]
        
        fusion_isoforms_indices = [[], []] 
        for read_idx in read_idx_ls:
            gname_idx = gene_names_ls.index(gname[read_idx])
            readname = fusion_isofroms_pair[read_idx][0]
            #start_isoform_idx = sum(num_isoforms_ls[:gname_idx])
            #end_isoform_idx = sum(num_isoforms_ls[:(gname_idx+1)])
            if (fusion_isoforms_map.has_key(readname)):
                fusion_isoforms_indices[read_idx] = [isoforms_names_ls.index(iso_name) for iso_name in fusion_isoforms_map[readname]]
            else:
                fusion_isoforms_indices[read_idx] = []
        num_fusion_isoforms = 0
        isoform_names_ls = fusion_MLE_input_dict[fusion_gene_name][1].split()
        isoform_len_ls = fusion_MLE_input_dict[fusion_gene_name][3].split()
        for isoform_idx_1 in fusion_isoforms_indices[0]:
            for isoform_idx_2 in fusion_isoforms_indices[1]:
                fusion_isoform_name = isoform_names_ls[isoform_idx_2][:-2] + '/' + isoform_names_ls[isoform_idx_1][:-2] # Remove .f from end
                if (fusion_isoform_name in processed_fusion_isoforms):  # Check both order
                    continue
                fusion_isoform_name = isoform_names_ls[isoform_idx_1][:-2] + '/' + isoform_names_ls[isoform_idx_2][:-2] # Remove .f from end
                if (fusion_isoform_name in processed_fusion_isoforms):
                    continue
                processed_fusion_isoforms.append(fusion_isoform_name)
                fusion_isoform_ls = [0] * num_gene_regions
                isoform_1 = [int(i) for i in fusion_MLE_input_dict[fusion_gene_name][7 + isoform_idx_1].split()]
                isoform_2 = [int(i) for i in fusion_MLE_input_dict[fusion_gene_name][7 + isoform_idx_2].split()]
                for i in range(len(isoform_1)):
                    fusion_isoform_ls[i] = max(isoform_1[i], isoform_2[i])
                for region in fusion_isofroms_pair[2]:
                    region_idx = gene_regions_ls.index(region)
                    fusion_isoform_ls[region_idx] = 1
                fusion_isoform_line = '\t'.join([str(i) for i in fusion_isoform_ls])
                fusion_MLE_input_dict[fusion_gene_name].insert(7 + num_isoforms + num_fusion_isoforms, fusion_isoform_line)
                # Add new isoform name
                isoform_names_ls.append(fusion_isoform_name)  
                isoform_len_ls.append(str(int(isoform_len_ls[isoform_idx_1]) + int(isoform_len_ls[isoform_idx_2])))
                num_fusion_isoforms += 1
        fusion_MLE_input_dict[fusion_gene_name][1] = '\t'.join(isoform_names_ls)
        fusion_MLE_input_dict[fusion_gene_name][3] = '\t'.join(isoform_len_ls)

        # Update number of isoforms item
        line_ls = fusion_MLE_input_dict[fusion_gene_name][0].split()
        line_ls[1] += ','+str(num_fusion_isoforms)
        fusion_MLE_input_dict[fusion_gene_name][0] = '\t'.join(line_ls)

    fusion_MLE_input_dict = finalize_fusion_genes(fusion_MLE_input_dict)
        
    return fusion_MLE_input_dict

###
########
def parse_sam_file(sam_filename, valid_pseudorefnames):
    
    print valid_pseudorefnames
    sam_filenames_dict = {}  # skip ref_000
    ref_id = 1
    for ref_id in valid_pseudorefnames:
        sam_filenames_dict[ref_id] = open(sam_filename + '.' + str(ref_id), 'w')
    

    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        if (line[0] == '@'):
            continue
        line_ls = line.split('\t')

        cigar_field = line_ls[5]
        if (len(set(cigar_field) - valid_cigar) > 0):
            continue

        ref_id = int(line_ls[2][4:])  # skip ref_
        if (ref_id in valid_pseudorefnames):
            sam_filenames_dict[ref_id].write(line)

    sam_file.close()
    for ref_id in valid_pseudorefnames:
        sam_filenames_dict[ref_id].close()
    return sam_filenames_dict

###
#########
def parse_cigar(line_ls):
    
    cigar_field = line_ls[5]
    cigar_list = re.split(r'(M|N|I|D)', cigar_field)
    # Note: this code is copied from parseSAM script!
    read_len_list = []
    seg_len = 0
    read_len = 0 
    M = 1
    for idx in range(len(cigar_list)/2):
        if (cigar_list[2 * idx + 1] == 'M'):
            if (M == 0):  # Mode is changed
                read_len_list.append(seg_len)
                seg_len = 0
            seg_len += int(cigar_list[2 * idx])
            read_len += int(cigar_list[2 * idx])
            M = 1
        elif (cigar_list[2 * idx + 1] == 'N'):
            if (M == 1):  # Mode is changed
                read_len_list.append(seg_len)
                seg_len = 0
            seg_len += int(cigar_list[2 * idx])
            M = 0
        elif (cigar_list[2 * idx + 1] == 'D'):  # Deletion from reference
            if (M == 0):  # Mode is changed
                read_len_list.append(seg_len)
                seg_len = 0
            seg_len += int(cigar_list[2 * idx])
        elif (cigar_list[2 * idx + 1] == 'I'):  # Insertion in reference
            if (M == 0):  # Mode is changed
                read_len_list.append(seg_len)
                seg_len = 0
            read_len +=  int(cigar_list[2 * idx])  
    read_len_list.append(seg_len)                
    
    if (abs(read_len - READ_LEN) > read_len_margin):
        read_len_list = []
    elif ((read_len_list[0] < READ_JUNC_MIN_MAP_LEN) or (read_len_list[-1] < READ_JUNC_MIN_MAP_LEN)):
        read_len_list = []
    
    return read_len_list
    
###
#########
def add_fusion_reads(fusion_MLE_input_dict, fusion_isofroms_pair_ls, fusion_gnames_map,
                     pseudo_ref_filename, sam_filename, compreadname2readname_mapping):

    
    num_aligned_fusion_reads = 0  # This would be added to header line of MLE_input.txt file
    pseudo_refname_dict =  {}  # ref_000 does not exits
    readid_2_refname_dict =  {}
    pseudo_ref_file = open(pseudo_ref_filename, 'r')
    ref_id = 1
    while (True):
        line = pseudo_ref_file.readline().strip()
        if not line:
            break
        pseudo_ref_file.readline()
        
        pseudo_refname_dict[line[1:]] = ref_id
        ref_id += 1
    pseudo_ref_file.close()

    readid_2_refname_dict =  {}
    pseudo_ref_file = open(pseudo_ref_filename + '.info', 'r')
    for line in pseudo_ref_file:
        fields = line.strip().split()
        readid_2_refname_dict[fields[0]] = fields[1:]
    pseudo_ref_file.close()
    
    valid_refnames = set()
    valid_pseudorefnames = set()
    for fusion_isofroms_pair in fusion_isofroms_pair_ls:
        if (len(fusion_isofroms_pair) < 2):
            continue # Not valid fusion gene
        # Note: because of removing redundancy, refnames might not match. That is the reason we check original rednames 
        origin_readnames = compreadname2readname_mapping[fusion_isofroms_pair[0][0]][fusion_isofroms_pair[1][0]]

        refname = set()
        for origin_readid in origin_readnames:
            refname.add(readid_2_refname_dict[origin_readid][0])
        if (len(refname) != 1):
            print "Err: Unexpected to have different refnames for fusion gene segments (expression level result could be invalid)" 
            print "\t" + str(fusion_isofroms_pair)
            print "\t" + str(refname)
        else:
            gname_1 = fusion_isofroms_pair[0][3]
            if (not fusion_gnames_map.has_key(gname_1)):
                continue
            valid_refnames.add(list(refname)[0])
            valid_pseudorefnames.add(pseudo_refname_dict[list(refname)[0]])

    sam_filenames_dict = parse_sam_file(sam_filename, valid_pseudorefnames)
    # fusion_isofroms_pair s could refer so same processed pair of genes, need to skip those
    processed_fusion_paires = []
    for fusion_isofroms_pair in fusion_isofroms_pair_ls:
        #print fusion_isofroms_pair
        if (len(fusion_isofroms_pair) < 2):
            continue # Not valid fusion gene

        [readname_1, fusion_1_pos] = re.split(r"\+|\-", (fusion_isofroms_pair[0][0].split("|"))[1])
        [readname_2, fusion_2_pos] = re.split(r"\+|\-", (fusion_isofroms_pair[1][0].split("|"))[1])
        gname_1 = fusion_isofroms_pair[0][3]
        gname_2 = fusion_isofroms_pair[1][3]
        if ((gname_1 + '/' + gname_2 in processed_fusion_paires) or
            (gname_2 + '/' + gname_1 in processed_fusion_paires)):
            pass
        else:
            processed_fusion_paires.append(gname_1 + '/' + gname_2)
        if (not fusion_gnames_map.has_key(gname_1)):
            continue
        if (not fusion_gnames_map.has_key(gname_2)):
            continue
        fusion_gene_name = fusion_gnames_map[gname_1]
        if (not fusion_MLE_input_dict.has_key(fusion_gene_name)):
            continue

        # Extract gene information
        gene_header_ls = fusion_MLE_input_dict[fusion_gene_name][0].split()
        gene_names_ls = gene_header_ls[0].split('/')
        chr_ls = gene_header_ls[2].split(',')
        point_names_ls = fusion_MLE_input_dict[fusion_gene_name][4].split()
        points_ls = [int(i) for i in fusion_MLE_input_dict[fusion_gene_name][5].split()]
        points_start_ls = [int(i) for i in gene_header_ls[3].split(',')]
        gene_regions_ls = fusion_MLE_input_dict[fusion_gene_name][6].split()
        regions_read_count_ls = fusion_MLE_input_dict[fusion_gene_name][-1].split()
        
        gname_1_idx = gene_names_ls.index(gname_1)
        gname_2_idx = gene_names_ls.index(gname_2)
        #Initialize the dict
        regions_count_dict = dict()
        for region_idx in range(int(gene_header_ls[4]), len(regions_read_count_ls)):
            regions_count_dict[gene_regions_ls[region_idx]] = 0
                
        #[LEFT, RIGHT] segment
        orig_id = list(compreadname2readname_mapping[fusion_isofroms_pair[0][0]][fusion_isofroms_pair[1][0]])[0] # Get one of them
        refname = readid_2_refname_dict[orig_id][0]
        #print pseudo_refname_dict[refname]
        if (pseudo_refname_dict[refname] not in valid_pseudorefnames):
            continue
        
        #if (pseudo_refname_dict[refname] != 10):
        #    """""""""DEBUGGING"""""""""
        #    continue
        #print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        #print fusion_isofroms_pair
        print fusion_isofroms_pair
        #print pseudo_refname_dict[refname]
        ref_name_ls = [refame_spl.split('_') for refame_spl in refname.split('/')]

        pseudo_ref_genome_points = [[int(ref_name_ls[0][1]), int(ref_name_ls[0][2])], [int(ref_name_ls[1][1]), int(ref_name_ls[1][2])]]
        pseudo_ref_offset = [0, pseudo_ref_genome_points[0][1] - pseudo_ref_genome_points[0][0]]
        pseudo_ref_dir = [ref_name_ls[0][3], ref_name_ls[1][3]]  

        # Need to detect first segment refers to which gene
        if ((ref_name_ls[0][0] == chr_ls[gname_1_idx]) and 
            (int(fusion_1_pos) > pseudo_ref_genome_points[0][0]) and
            (int(fusion_1_pos) < pseudo_ref_genome_points[0][1])):
            gname_idx_ls = [gname_1_idx, gname_2_idx]
            fusion_idx_ls = [int(fusion_1_pos), int(fusion_2_pos)] 
        else:
            gname_idx_ls = [gname_2_idx, gname_1_idx]
            fusion_idx_ls = [int(fusion_2_pos), int(fusion_1_pos)]
        sam_file = open(sam_filename + '.' + str(pseudo_refname_dict[refname]),'r')
        for read_line in sam_file:
            read_line_ls = read_line.split()
            #if (read_line_ls[0] != "SRR1107833.2819532"):
            #    """""""""DEBUGGING"""""""""
            #    continue
            #print read_line
            read_start_pos = int(read_line_ls[3])
            read_len_ls = parse_cigar(read_line_ls)
            if (len(read_len_ls) == 0) :
                continue
            region_str_ls = [[], []]
            read_len_idx = 0
            read_genome_pos = read_start_pos
            read_genome_pos_offset = 0  # offset for second segment
            valid_read = True

            for seg_idx in [0, 1]:
                #print region_str_ls
                if (seg_idx == 1):
                    if (read_len_idx > len(read_len_ls)):
                        valid_read = False
                        break
                    #print read_len_idx
                    #print read_len_ls
                    #print read_start_pos + sum(read_len_ls[:read_len_idx])
                    #print pseudo_ref_offset[1]
                    # For second segment we should have proceeded to second part of reference seq
                    if ((read_start_pos + sum(read_len_ls[:read_len_idx])) < pseudo_ref_offset[1]):
                        valid_read = False
                        break
                    read_genome_pos = read_start_pos + sum(read_len_ls[:read_len_idx]) - pseudo_ref_offset[1]
                if (valid_read == False):
                    break
                if (pseudo_ref_dir[seg_idx] == "+"):
                    read_genome_pos += pseudo_ref_genome_points[seg_idx][0]
                    point_idx_ls = range(points_start_ls[gname_idx_ls[seg_idx]], points_start_ls[gname_idx_ls[seg_idx]+1])  # List of point indices for this gene
                    # Note: the resulter region might not be valid which in case it wont be mapped to any gene region
                    for point_idx in point_idx_ls:
                        #print points_ls[point_idx], read_genome_pos
                        if (points_ls[point_idx] > pseudo_ref_genome_points[seg_idx][1]):
                            break

                        if (points_ls[point_idx] < read_genome_pos):
                            continue
                        elif (points_ls[point_idx] == read_genome_pos):
                            region_str_ls[seg_idx].append(point_names_ls[point_idx])
                        else:
                            if ((read_genome_pos + read_len_ls[read_len_idx] - 1) > points_ls[point_idx]):
                                region_str_ls[seg_idx].append(point_names_ls[point_idx])
                            elif ((read_genome_pos + read_len_ls[read_len_idx] - 1) == points_ls[point_idx]):
                                region_str_ls[seg_idx].append(point_names_ls[point_idx])
                                read_genome_pos += read_len_ls[read_len_idx] - 1  #do not include gap point
                                read_len_idx += 2
                                if (read_len_idx > len(read_len_ls)):
                                    break
                                # Add previous gap
                                read_genome_pos += read_len_ls[read_len_idx - 1] + 1
                                if (seg_idx == 0):
                                    if (read_genome_pos > fusion_idx_ls[0]):
                                        break
                            else:
                                if (seg_idx == 1):
                                    # Last block partialy extened within exon block
                                    if ((read_len_idx + 1) == len(read_len_ls)):
                                        if ((int(point_names_ls[point_idx][1:]) % 2) == 1):
                                            break
                                        elif ((point_idx > points_start_ls[gname_idx_ls[seg_idx]]) and
                                              ((points_ls[point_idx] - points_ls[point_idx-1]) == 1)):
                                            break
                                        else:
                                            valid_read = False
                                            break
                                else:
                                    valid_read = False
                                    break
                else:
                    read_genome_pos = pseudo_ref_genome_points[seg_idx][1] + 1 - read_genome_pos
                    point_idx_ls = list(reversed(range(points_start_ls[gname_idx_ls[seg_idx]], points_start_ls[gname_idx_ls[seg_idx]+1])))
                    #print point_idx_ls
                    for point_idx in point_idx_ls:
                        #print points_ls[point_idx], read_genome_pos
                        if (points_ls[point_idx] < pseudo_ref_genome_points[seg_idx][0]):
                            break
                        if (points_ls[point_idx] > read_genome_pos):
                            continue
                        elif (points_ls[point_idx] == read_genome_pos):
                            region_str_ls[seg_idx].append(point_names_ls[point_idx])
                        else:
                            if ((read_genome_pos - (read_len_ls[read_len_idx] - 1)) < points_ls[point_idx]):
                                region_str_ls[seg_idx].append(point_names_ls[point_idx])
                            elif ((read_genome_pos - (read_len_ls[read_len_idx] - 1)) == points_ls[point_idx]):
                                region_str_ls[seg_idx].append(point_names_ls[point_idx])
                                read_genome_pos -= read_len_ls[read_len_idx] - 1
                                read_len_idx += 2
                                if (read_len_idx > len(read_len_ls)):
                                    break
                                read_genome_pos -=  read_len_ls[read_len_idx - 1] + 1 # Compensate for above -1
                                if (seg_idx == 0):
                                    if (read_genome_pos < fusion_idx_ls[0]):
                                        break
                            else:
                                if (seg_idx == 1):
                                    # Last seg partialy overlap with exon
                                    if ((read_len_idx + 1) == len(read_len_ls)):
                                        if ((int(point_names_ls[point_idx][1:]) % 2) == 0):
                                            break
                                        elif (((point_idx + 1) < points_start_ls[gname_idx_ls[seg_idx] + 1]) and
                                              ((points_ls[point_idx+1] - points_ls[point_idx]) == 1)):
                                            break
                                    
                                        else:
                                            valid_read = False
                                            break
                                valid_read = False
                                break
            ##print valid_read

            if (valid_read == False):
                continue
            region_str_ls = region_str_ls[0] + region_str_ls[1]
            if (gname_idx_ls[0] < gname_idx_ls[1]):
                region_str = '-'.join(region_str_ls)
            else:
                region_str = '-'.join(list(reversed(region_str_ls)))

            if (regions_count_dict.has_key(region_str)):
                regions_count_dict[region_str] += 1

        #Update the read-count
        #print regions_count_dict
        for region_idx in range(int(gene_header_ls[4]), len(regions_read_count_ls)):
            if (int(regions_read_count_ls[region_idx]) == 0):  # To avoid replicates, it should be parse only once
                regions_read_count_ls[region_idx] = regions_count_dict[gene_regions_ls[region_idx]]
                num_aligned_fusion_reads += regions_count_dict[gene_regions_ls[region_idx]]
        fusion_MLE_input_dict[fusion_gene_name][-1] = '\t'.join([str(i) for i in regions_read_count_ls])
        print regions_count_dict
    return [fusion_MLE_input_dict, num_aligned_fusion_reads] 
            
            
def print_MLE_input_file(fusion_MLE_input_dict, fusion_gnames_map, num_aligned_fusion_reads,
                         input_filename, output_filename):


    input_file = open(input_filename, 'r')
    output_file = open(output_filename, 'w')
    output_file_ = open(output_filename+'_', 'w')  # Keep a seperate copy of fusion isoforms
    fusion_gene_names_ls = fusion_MLE_input_dict.keys()
    line_index = 0;
    num_aligned_reads = int(input_file.readline()) + num_aligned_fusion_reads 
    output_file.write(str(num_aligned_reads) + '\n')
    wr_flag = True
    print fusion_gene_names_ls
    for line in input_file:
        fields = line.split()
        if (line_index == 0):
            gname = fields[0]
            num_isoforms = int(fields[1])
            rname = fields[2]
            if (fusion_gnames_map.has_key(gname)):
                fusion_gname = fusion_gnames_map[gname]
                if (fusion_gname in fusion_gene_names_ls):
                    wr_flag = False
        line_index += 1
        if (wr_flag):
            output_file.write(line)
        if (line_index == (9 + num_isoforms)):
            line_index = 0
            wr_flag = True
    
    for fusion_gname in fusion_gene_names_ls:
        for fusion_line in fusion_MLE_input_dict[fusion_gname]:
            output_file.write(fusion_line.strip() + '\n')
            output_file_.write(fusion_line.strip() + '\n')

    input_file.close()
    output_file.close()
    output_file_.close()
    


### Main
##########
def main():
    
    fusion_compatible_filename = sys.argv[1]
    fusion_gpd_filename = sys.argv[2]
    readmapping_filename = sys.argv[3]
    isoforms_readnames_filename = sys.argv[4]
    input_filename = sys.argv[5]
    
    sam_filename = sys.argv[6]
    pseudo_ref_filename = sys.argv[7]
    output_filename = sys.argv[8]
    
    global READ_LEN, READ_JUNC_MIN_MAP_LEN
    READ_LEN = int(sys.argv[9])
    READ_JUNC_MIN_MAP_LEN = int(sys.argv[10])
    
    [fusion_isofroms_pair_ls, MLE_input_dict, compreadname2readname_mapping] = parse_fusion_genenames(fusion_compatible_filename,
                                                                                                      fusion_gpd_filename, readmapping_filename)
    MLE_input_dict = get_fusion_genomes(input_filename, MLE_input_dict)
    
    #fusion_isofroms_pair_ls = map_id2genome(readmapping_filename, fusion_compatible_filename,
    #                               fusion_isofroms_pair_ls, MLE_input_dict)
    print ">>>fusion_isofroms_pair_ls"
    for i in fusion_isofroms_pair_ls:
        print i   

    [fusion_MLE_input_dict, fusion_gnames_map] = merge_fusion_genes(MLE_input_dict, fusion_isofroms_pair_ls)
    print ">>>fusion_gnames_map"
    for key in fusion_gnames_map.keys():
        print key +": " + fusion_gnames_map[key]
    
        
    fusion_MLE_input_dict = add_fusion_regions(fusion_MLE_input_dict, fusion_gnames_map, fusion_isofroms_pair_ls)
    fusion_MLE_input_dict = add_fusion_isoforms(fusion_MLE_input_dict, fusion_gnames_map, fusion_isofroms_pair_ls,
                                                isoforms_readnames_filename)

    [fusion_MLE_input_dict, num_aligned_fusion_reads] = add_fusion_reads(fusion_MLE_input_dict, fusion_isofroms_pair_ls, fusion_gnames_map,
                                                                         pseudo_ref_filename, sam_filename, compreadname2readname_mapping)

    print_MLE_input_file(fusion_MLE_input_dict, fusion_gnames_map, num_aligned_fusion_reads,
                         input_filename, output_filename)
    """
    for i in compreadname2readname_mapping.keys():
        for j in compreadname2readname_mapping[i].keys():
            print i + " " + j
            print compreadname2readname_mapping[i][j]
    """
    """
    print ">>>Original genes"
    for i in MLE_input_dict.keys():
        for j in MLE_input_dict[i]:
            print j
        print "-----"

    print ">>>Fusion genes"
    for i in fusion_MLE_input_dict.keys():
        for j in fusion_MLE_input_dict[i]:
            print j
        print "-----"
    """
if __name__ == '__main__':
    main()





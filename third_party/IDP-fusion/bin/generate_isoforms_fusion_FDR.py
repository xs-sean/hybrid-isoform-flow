#!/usr/bin/python
import sys
import re
import string
#from constants import *
import re



    
### Store fusion gene information from MLE input file
#########
def get_fusion_genomes(input_filename, isoforms_name_filename):
    
    MLE_input_dict = dict()
    isoforms_names_ls = []
    isoforms_single_names_ls = []
    isoforms_name_file = open(isoforms_name_filename, 'r')
    for line in isoforms_name_file:
        if ("/" in line):
            isoforms_single_names_ls += line.strip().split("/")
            isoforms_names_ls.append(line.strip())
    isoforms_name_file.close()
    num_reads = 0
    
    input_file = open(input_filename, 'r')
    line_index = 0;
    wr_flag = False
    for line in input_file:
        fields = line.split()
        if (line_index == 0):
            gname = fields[0]
            num_isoforms = int(fields[1])
            rname = fields[2]
            prev_line = line
        elif (line_index == 1):
            for isoform in line.strip().split():
                if (isoform in isoforms_names_ls):
                     MLE_input_dict[gname] = []
                     MLE_input_dict[gname].append(prev_line)
                     wr_flag = True
                     break
        line_index += 1
        if (wr_flag):
            MLE_input_dict[gname].append(line.strip())
        if (line_index == (9 + num_isoforms)):
            line_index = 0
            wr_flag = False
            
    input_file.close()

    # Keep region with supporting isoforms
    gnames = MLE_input_dict.keys()
    for gname in gnames:

        isoforms_ls = MLE_input_dict[gname][1].strip().split()
        isoforms_flag_ls = [(i in isoforms_names_ls) for i in isoforms_ls]
        num_isoforms = len(isoforms_ls)
        num_regions = len(MLE_input_dict[gname][6].strip().split())
        regions_flag_ls = [0] * num_regions
        for isoform_idx in range(num_isoforms):
            if (isoforms_flag_ls[isoform_idx]):
                line_ls = [int(i) for i in MLE_input_dict[gname][7 + isoform_idx].strip().split()]
                regions_flag_ls = [regions_flag_ls[i] + line_ls[i] for i in range(num_regions)]
            
        for line_idx in [6, 8 + isoform_idx, 9 + isoform_idx]:
            line_ls = MLE_input_dict[gname][line_idx].strip().split()
            new_line_ls = [line_ls[i] for i in range(num_regions) if regions_flag_ls[i] > 0]
            MLE_input_dict[gname][line_idx] = "\t".join(new_line_ls)

        # Update region columns
        
        for isoform_idx in range(num_isoforms-1, -1, -1):
            if (isoforms_flag_ls[isoform_idx]):
                line_ls = MLE_input_dict[gname][7 + isoform_idx].strip().split()
                new_line_ls = [line_ls[i] for i in range(num_regions) if regions_flag_ls[i] > 0]
                MLE_input_dict[gname][7 + isoform_idx] = "\t".join(new_line_ls)
            else:
                del MLE_input_dict[gname][7 + isoform_idx]
            
        # Update list of isoforms
        isoforms_ls = MLE_input_dict[gname][1].strip().split()
        for line_idx in [1, 2, 3]:
            line_ls = MLE_input_dict[gname][line_idx].strip().split()
            new_line_ls = [line_ls[i] for i in range(num_isoforms) if isoforms_ls[i] in isoforms_names_ls]
            MLE_input_dict[gname][line_idx] = "\t".join(new_line_ls)
        
        # Update number of isoforms
        line_ls = MLE_input_dict[gname][0].strip().split()
        line_ls[1] = str(sum(isoforms_flag_ls))
        MLE_input_dict[gname][0] = "\t".join(line_ls)
        
        num_reads_ls = [int(i) for i in  MLE_input_dict[gname][-1].strip().split()]
        num_reads += sum(num_reads_ls)


    return [MLE_input_dict, isoforms_single_names_ls, num_reads]



def print_MLE_input_file(input_filename, output_filename,
                         isoforms_names_filename, MLE_input_dict, num_reads):


    input_file = open(input_filename, 'r')
    output_file = open(output_filename, 'w')
    
    num_reads_origin = int(input_file.readline())
    output_file.write(str(num_reads_origin + num_reads) + "\n")
    for line in input_file:
        output_file.write(line)
        
    fusion_gene_names_ls = MLE_input_dict.keys()
    for fusion_gname in fusion_gene_names_ls:
        for fusion_line in MLE_input_dict[fusion_gname]:
            output_file.write(fusion_line.strip() + '\n')
            
    input_file.close()
    output_file.close()
    

def print_predicted_gpd_file(input_gpd_filename, output_gpd_filename, isoforms_names_ls):
    
    input_gpd_file = open(input_gpd_filename, 'r')
    output_gpd_file = open(output_gpd_filename, 'w')
    print isoforms_names_ls
    for line in input_gpd_file:
        fields = line.split()
        if ((fields[1][:-2]) in isoforms_names_ls):
            output_gpd_file.write(line)
    
    input_gpd_file.close()
    output_gpd_file.close()
    
    

### Main
##########
def main():
    
    fusion_parseSAM_filename = sys.argv[1]
    selected_fusion_isoforms_filename = sys.argv[2]
    input_parseSAM_filename = sys.argv[3]
    output_parseSAM_filename = sys.argv[4]
    input_gpd_filename = sys.argv[5]
    output_gpd_filename = sys.argv[6]   
     
    [MLE_input_dict, isoforms_names_ls, num_reads] = get_fusion_genomes(fusion_parseSAM_filename, selected_fusion_isoforms_filename)
    

    print_MLE_input_file(input_parseSAM_filename, output_parseSAM_filename, selected_fusion_isoforms_filename, 
                         MLE_input_dict, num_reads)
    print_predicted_gpd_file(input_gpd_filename, output_gpd_filename, isoforms_names_ls)
    
if __name__ == '__main__':
    main()





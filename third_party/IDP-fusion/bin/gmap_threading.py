#!/usr/bin/python
import sys
import struct
import os
import string

################################################################################

if len(sys.argv) >= 7:
    gmap_path = sys.argv[1]
    gmap_option = ' '.join(sys.argv[2:-3])
    LR_pathfilename = sys.argv[-3]
    gmap_index_pathfoldername = sys.argv[-2]
    output_pathfilename = sys.argv[-1]
    
else:
    print("usage: python2.6 gmap_threading.py -f 1 -D ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/ -d genome/ /usr/bin/python intact_SM.fa intact_SM.fa.psl")
    print("or ./blat_threading.py p -t=DNA -q=DNA ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /usr/bin/python intact_SM.fa intact_SM.fa.psl")
    sys.exit(1)
################################################################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    if path == "/":
        path = "./"
    return path, filename
################################################################################
LR_path, LR_filename = GetPathAndName(LR_pathfilename)
output_path, output_filename = GetPathAndName(output_pathfilename)
gmap_folder, gmap_index = GetPathAndName(gmap_index_pathfoldername)

################################################################################
#print "Warning: gmap threading is disabled"
#exit(1)
LR = open(LR_pathfilename,'r')
LR_NR = 0
for line in LR:
    LR_NR+=1
LR.close()

##########################################

gmap_LR_cmd = gmap_path + " " + gmap_option + ' ' + " -D " + gmap_folder + " -d " + gmap_index + " " + LR_pathfilename + ' > ' + output_pathfilename  
print gmap_LR_cmd
os.system(gmap_LR_cmd)
	
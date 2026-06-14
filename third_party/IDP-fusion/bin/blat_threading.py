#!/usr/bin/python
import sys
import struct
import os
import threading
import string

################################################################################

if len(sys.argv) >= 7: 
    blat_path = sys.argv[1]
    Nthread1 = int(sys.argv[2])
    blat_option = ' '.join(sys.argv[3:-2])
    LR_pathfilename = sys.argv[-2]
    output_pathfilename = sys.argv[-1]
    
else:
    print("usage: python blat_threading.py p -t=DNA -q=DNA ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa intact_SM.fa intact_SM.fa.psl")
    print("or ./blat_threading.py p -t=DNA -q=DNA ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa intact_SM.fa intact_SM.fa.psl")
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
################################################################################
#print "Warning: blat threading is disabled"
#exit(1)
LR = open(LR_pathfilename,'r')
LR_NR = 0
for line in LR:
    LR_NR+=1
LR.close()

Nsplitline = 1 + (LR_NR/Nthread1)
if Nsplitline%2==1:
    Nsplitline +=1
ext_ls=[]
j=0
k=0
i=0

Nthread = LR_NR/Nsplitline
if (LR_NR%Nsplitline > 0):
    Nthread += 1
Nthread1 = min(Nthread1, Nthread)
while i < Nthread1:
    ext_ls.append( '.' + string.lowercase[j] + string.lowercase[k] )
    k+=1
    if k==26:
        j+=1
        k=0
    i+=1

print "===split LR:==="    
splitLR_cmd = "split -l " + str(Nsplitline) + " " + LR_pathfilename + " " + output_path + LR_filename +"."
print splitLR_cmd
os.system(splitLR_cmd)

##########################################
print "===Run blat:==="    

i=0
T_blat_LR_ls = []
for ext in ext_ls:
    blat_LR_cmd = blat_path + " " + blat_option + ' ' + output_path + LR_filename + ext + ' ' + output_path + LR_filename + ext + ".psl"
    print blat_LR_cmd
    T_blat_LR_ls.append( threading.Thread(target=os.system, args=(blat_LR_cmd,)) )
    T_blat_LR_ls[i].start()
    i+=1

for T in T_blat_LR_ls:
    T.join()

cat_psl_cmd = "cat "
rm_LR_cmd = "rm "
rm_LRpsl_cmd = "rm "
for ext in ext_ls:
    cat_psl_cmd = cat_psl_cmd + output_path + LR_filename + ext + ".psl "
    rm_LR_cmd = rm_LR_cmd + output_path + LR_filename + ext + ' '
    rm_LRpsl_cmd = rm_LRpsl_cmd + output_path + LR_filename + ext + ".psl "
cat_psl_cmd = cat_psl_cmd + " > " + output_pathfilename    
print cat_psl_cmd
print rm_LR_cmd
print rm_LRpsl_cmd

os.system(cat_psl_cmd)
os.system(rm_LR_cmd)
os.system(rm_LRpsl_cmd)

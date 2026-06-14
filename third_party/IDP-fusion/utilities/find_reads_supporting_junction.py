#!/usr/bin/python
import sys, os, inspect
import re

#pre: 1.) the temporary directory created by an IDP-fusion run, 
#     2.) fusion point coordinates in the format chr1:1102299+/chr2:220111-
#         where the coordinate on either side is the 1-indexed coordiante of the
#         last base in the transcript before the fusion begins
#         '+' on the left side of '/' means that the junction is on the right side 
#         '-' means it is on the left 
#         '+' on the right of '/' means that the junction on that side is on the left
#         '-' means its on the right
#     3.) Range .. how far to look from the junction for supporting long reads
#     4.) output directory .. where to write outputs to
# Post:  Creates files in the output directory
#   fusion_coordiante.txt contains the fusion coordiante that was input followed by the a tab then the range to search for long reads
#   long_read_query_locations.txt says where the left and right sides fall on
#        the query sequences
#   long_read_left.bed say where the left side of junctions long reads map to the genome
#   long_read_right.bed say where the right side of the junctions long reads map tot he genome
#   short_read_left.bed say where the left side of the junctions short reads map
#   short_read_right.bed say where the right side of the junctions short reads map to the genome

#bring in the folder to the path for our modules
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],"../bin")))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import psl_basics
import genepred_basics

def main():

  if len(sys.argv) != 5:
    print sys.argv[0] + ' <IDP temp dir> <junction (i.e. chr1:1000-/chr2:1000+)> <range> <output dir>'
    sys.exit()
  inputdir = sys.argv[1].rstrip('/')
  psl_filename = inputdir+'/'+'LR.psl_fusion_pair_filtered'
  map_filename = inputdir+'/uniqueness/fusion_read.map'
  fasta_filename = inputdir+'/uniqueness/unmapped_shortread.fa'
  txn_filename = inputdir+'/uniqueness/txn.map'
  junction_abbreviation = sys.argv[2]
  myrange = int(sys.argv[3])
  outdir = sys.argv[4].rstrip('/')
  if not os.path.isdir(outdir):
    print "make directory "+outdir
    os.mkdir(outdir)
  of1 = open(outdir+"/fusion_coordiante.txt",'w')
  of1.write(junction_abbreviation+"\t"+str(myrange)+"\n")
  of1.close()
  m = re.match('([^:]+):(\d+)([+-])\/([^:]+):(\d+)([+-])',junction_abbreviation)
  chr1 = m.group(1)
  coo1 = int(m.group(2))
  dir1 = m.group(3)
  chr2 = m.group(4)
  coo2 = int(m.group(5))
  dir2 = m.group(6)
  of = open(outdir+"/long_read_query_locations.txt",'w')
  oleft = open(outdir+"/long_read_left.bed",'w')
  oright = open(outdir+"/long_read_right.bed",'w')
  oleftgpd = open(outdir+"/long_read_left.gpd",'w')
  orightgpd = open(outdir+"/long_read_right.gpd",'w')
  # Work through psl file for long reads
  lrcnt = 0
  with open(psl_filename) as f:
    while True:
      l1 = f.readline().rstrip()
      if not l1:
        break
      l2 = f.readline().rstrip()
      if not l2:
        break
      e1 = psl_basics.read_psl_entry(l1)
      e2 = psl_basics.read_psl_entry(l2)
      g1 = genepred_basics.smooth_gaps(genepred_basics.genepred_line_to_dictionary(psl_basics.convert_entry_to_genepred_line(e1)),30)
      g2 = genepred_basics.smooth_gaps(genepred_basics.genepred_line_to_dictionary(psl_basics.convert_entry_to_genepred_line(e2)),30)
      if check_coordiantes(e1,e2,chr1,coo1,dir1,chr2,coo2,dir2,myrange):
        oleftgpd.write(genepred_basics.genepred_entry_to_genepred_line(g1)+"\n")
        orightgpd.write(genepred_basics.genepred_entry_to_genepred_line(g2)+"\n")
        lrcnt += 1
        of.write(e1["qName"] + "\t" + str(e1['qStart']+1) + "\t" + \
                 str(e1['qEnd']) + "\t" + dir1 + "\t" + \
                 str(e2['qStart']+1) + "\t" + str(e2['qEnd']) + "\t" + dir2 + "\n")
        for i in range(0,g1['exonCount']):
          if g1["exonEnds"][i]-g1["exonStarts"][i] >= 30:
            oleft.write(g1["chrom"] + "\t" + str(g1["exonStarts"][i]+1) + "\t" + \
                        str(g1["exonEnds"][i]) + "\t" + e1["qName"] + "\n")
        for i in range(0,g2['exonCount']):
          if g2["exonEnds"][i]-g2["exonStarts"][i] >= 30:
            oright.write(g2["chrom"] + "\t" + str(g2["exonStarts"][i]+1) + "\t" + \
                        str(g2["exonEnds"][i]) + "\t" + e1["qName"] + "\n")
      elif check_coordiantes(e1,e2,chr2,coo2,opposite(dir2),chr1,coo1,opposite(dir1),myrange):
        oleftgpd.write(genepred_basics.genepred_entry_to_genepred_line(g2)+"\n")
        orightgpd.write(genepred_basics.genepred_entry_to_genepred_line(g1)+"\n")
        lrcnt += 1
        of.write(e1["qName"] + "\t" + str(e1['qStart']+1) + "\t" + \
                 str(e1['qEnd']) + "\t" + opposite(dir2) + "\t" + \
                 str(e2['qStart']+1) + "\t" + str(e2['qEnd']) + "\t" + opposite(dir1) + "\n")
        for i in range(0,g1['exonCount']):
          if g1["exonEnds"][i]-g1["exonStarts"][i] >= 30:
            oright.write(g1["chrom"] + "\t" + str(g1["exonStarts"][i]+1) + "\t" + \
                         str(g1["exonEnds"][i]) + "\t" + e1["qName"] + "\n")
        for i in range(0,g2['exonCount']):
          if g2["exonEnds"][i]-g2["exonStarts"][i] >= 30:
            oleft.write(g2["chrom"] + "\t" + str(g2["exonStarts"][i]+1) + "\t" + \
                        str(g2["exonEnds"][i]) + "\t" + e1["qName"] + "\n")
  of.close()
  oleft.close()
  oright.close()
  oleftgpd.close()
  orightgpd.close()
  print str(lrcnt) + " long reads found supporting the fusion"
  #Work through fusion read map for short reads
  rnames = {}
  seenhit = {}
  with open(map_filename) as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      srname = f[0]
      loc = f[2]
      m = re.search('^([^:]+):.*\/([^:]+):',loc)
      srchr1 = m.group(1)
      srchr2 = m.group(2)
      m = re.search('[,:](-?\d+)-(-?\d+)([+-])\/[^:]+:(-?\d+)-(-?\d+),?.*([+-])',loc)
      srcoo1start = int(m.group(1))
      srcoo1finish = int(m.group(2))
      srdir1 = m.group(3)
      srcoo2start = int(m.group(4))
      srcoo2finish = int(m.group(5))
      srdir2 = m.group(6)
      m = re.search
      srcooleft = srcoo1finish
      if srdir1 == '-':
        srcooleft = srcoo1start
      srcooright = srcoo2start
      if srdir2 == '-':
        srcooright = srcoo2finish
      #print srchr1 + "\t" + srchr2 + "\t" + str(srcooleft) + "\t" + str(srcooright)
      if srdir1 == dir1 and srchr1 == chr1 and srdir2 == dir2 and srchr2 == chr2 and srcooleft == coo1 and srcooright == coo2: 
        rnames[srname] = {}
        rnames[srname]['left'] = srchr1 + "\t" + str(srcoo1start) + "\t" + str(srcoo1finish) + "\t" + srname
        rnames[srname]['right'] = srchr2 + "\t" + str(srcoo2start) + "\t" + str(srcoo2finish) + "\t" + srname
        if srname not in seenhit:
          seenhit[srname] = 0
        seenhit[srname] += 1
      if srdir1 == opposite(dir1) and srchr1 == chr2 and srdir2 == opposite(dir2) and srchr2 == chr1 and srcooleft == coo2 and srcooright == coo1: 
        rnames[srname] = {}
        rnames[srname]['left'] = srchr2 + "\t" + str(srcoo2start) + "\t" + str(srcoo2finish) + "\t" + srname
        rnames[srname]['right'] = srchr1 + "\t" + str(srcoo1start) + "\t" + str(srcoo1finish) + "\t" + srname
        if srname not in seenhit:
          seenhit[srname] = 0
        seenhit[srname] += 1
  print "found "+str(len(rnames))+" short reads"
  for srname in seenhit:
    if seenhit[srname] > 1:
      print "removing " + srname
      del rnames[srname]
  print "found "+str(len(rnames))+" short reads with no multihits among fusions"
  validreads = {}
  with open(fasta_filename) as inf:
    for line in inf:
      m = re.match('^>(.*)$',line.rstrip())
      if m:
        validreads[m.group(1)] = 1
  namelist = rnames.keys()
  for rname in namelist:
    if rname not in validreads:
      print "removing " + rname
      del rnames[rname]
  print "found "+str(len(rnames))+" unique short reads with no hits in the genome"
  with open(txn_filename) as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      if f[0] in rnames:
        print "removing " + f[0]
        del rnames[f[0]]
  print "found "+str(len(rnames))+" unique short reads with no hits in the genome or transcriptome"
  oleft = open(outdir+"/short_read_left.bed",'w')
  oright = open(outdir+"/short_read_right.bed",'w')
  for rname in rnames:
    oleft.write(rnames[rname]['left']+"\n")
    oright.write(rnames[rname]['right']+"\n")
  oleft.close()
  oright.close()

def near(c1,c2,myrange):
  if abs(c1-c2) <= myrange:
    return True
  return False

def check_coordiantes(e1,e2,chr1,coo1,dir1,chr2,coo2,dir2,myrange):
  if e1['tName'] == chr1 and e2['tName'] == chr2:
    if (dir1 == '+' and \
      coo1 <= e1['tEnd']+myrange and \
      coo1 >= e1['tEnd']-myrange) or \
      (dir1 == '-' and \
      coo1 <= e1['tStart']+myrange and \
      coo1 >= e1['tStart']-myrange) and \
      (dir2 == '+' and \
      coo2 <= e2['tStart']+myrange and \
      coo2 >= e2['tStart']-myrange) or \
      (dir2 == '-'	and \
      coo2 <= e2['tEnd']+myrange and \
      coo2 >= e2['tEnd']-myrange):
        return True
  return False

def opposite(sign):
  if sign=='+':
    return '-'
  return '+'

main()

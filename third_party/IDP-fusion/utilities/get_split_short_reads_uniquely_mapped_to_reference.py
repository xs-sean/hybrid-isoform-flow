#!/usr/bin/python
import sys, os, inspect
import re


#bring in the folder to the path for our modules
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],"../bin")))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import sam_basics

def main():

  if len(sys.argv) != 2:
    print sys.argv[0] + ' <IDP temp dir>'
    sys.exit()
  inputdir = sys.argv[1].rstrip('/')
  sam_filename = inputdir+'/uniqueness/genome.sam'
  read_filename = inputdir+'/reference_read_counts.txt'
  nonuniquereads = set()
  i = 0
  with open(read_filename) as inf:
    for line in inf:
      i+=1
      if i%1000000==0:
        print "read "+str(i)+ " reads"
      f = line.rstrip().split("\t")
      srname = f[0]
      cnt = int(f[1])
      if cnt != 1:
        nonuniquereads.add(srname)

  with open(sam_filename) as inf:
    for line in inf:
      line = line.rstrip()
      if sam_basics.is_header(line):
        print line
        continue
      d = sam_basics.sam_line_to_dictionary(line)
      if d['qname'] in nonuniquereads:
        continue
      print line

main()

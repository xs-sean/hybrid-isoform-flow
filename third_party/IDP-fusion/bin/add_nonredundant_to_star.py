#!/usr/bin/python

import sys;
import re;

minimum_overlap = 12

if len(sys.argv) == 4:
  star_sam_filename = sys.argv[1]
  junctions_filename = sys.argv[2]
  minimum_overlap = int(sys.argv[3])
else:
  print("usage:add_nonredundant_to_star.py Aligned.out.sam SJ.out.tab min-overlap-len")
  sys.exit(1)
j = {}
# read the junction file and coordinates into a data structure
# indexed by chromosome
with open(junctions_filename) as infile:
  for line in infile:
    f = line.rstrip().split("\t")
    if not f[0] in j:
      j[f[0]] = []
    t = {}
    t['first'] = f[1]
    t['last'] = f[2]
    t['nrcount'] = 0 # this is what we are calculating.  nonredundant count.
    t['fields']  = f
    j[f[0]].append(t)

# load the reads into nonredundant entries with a dictionary
r = {}
with open(star_sam_filename) as infile:
  for line in infile:
    f = line.rstrip().split("\t")
    a = re.compile("^([0-9]+)M([0-9]+)N([0-9]+)M") #get ready to skip if its an exact match
    if len(f) == 15 and a.match(f[5]) and f[2] in j:
      chrom = f[2]
      pos = f[3]
      # see if the sam line matches a dictionary entry
      m = a.match(f[5])
      leftcount = m.group(1)
      spacecount = m.group(2)
      rightcount = m.group(3)
      #only add it as a junction if it meets minimum overlap requirements
      if int(leftcount) >= minimum_overlap and int(rightcount) >= minimum_overlap:
        if not chrom in r:
          r[chrom] = {}
        v = "\t".join([pos, leftcount, spacecount, rightcount])
        r[chrom][v] = 1

for chrom in j:
  if chrom in r:
    for junc in j[chrom]:
      first = int(junc['first'])
      last = int(junc['last'])
      #print "\t".join(junc['fields'])
      for sam in r[chrom]:
        f = sam.split("\t")
        intron1 = int(f[0]) + int(f[1]) 
        intron2 = int(f[0]) + int(f[1]) + int(f[2]) - 1
        if first == intron1 and last == intron2:
          junc['nrcount'] += 1
     
for chrom in j:
  for entry in j[chrom]:
    print "\t".join(entry['fields']) + "\t" + str(entry['nrcount'])

#!/usr/bin/python
"""Parse cluster list.
python parse_member_list_may26.py > ../may26.k235.avglink.5ct.extant.collapsed.hclust.list.txt
"""
FNAME = "../may26.k235.avglink.5ct.extant.collapsed.hclust.csv"
clust = {}
for line in open(FNAME):
  g,c = line.strip().split(',')
  clust.setdefault(int(c),set()).add(g.strip('"'))
for c in sorted(clust):
  print "%d,%s" % (c, ",".join(sorted(clust[c])))
  

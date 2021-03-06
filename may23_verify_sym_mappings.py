#!/usr/bin/python
"""
good: 13632 remapped 187 not found 396 multiple names 3

Create a map from original gene symbol to HUGO remapped symbol if possible.
"""
import hugo_gene_symbols
H = hugo_gene_symbols.load()
FNAME = "E7307.original.row.names.txt"
FNAME_OUT = "E7307.genesym.remap.txt"

good = set()
other = {}
bad = set()
dupe = {}


fp_out = open(FNAME_OUT,"w")
for line in open(FNAME):
  s = line.strip()
  z = H.find_sym(s, allow_dupe=True)
  if s != z:
    if isinstance(z,set):
      dupe[s] = z
      s_map = s
    elif z:
      other[s] = z
      s_map = z
    else:
      bad.add(s)
      s_map = s
  else:
    good.add(s)
    s_map = s
  fp_out.write("%s\t%s\n" % (s,s_map))
fp_out.close()

print "good:", len(good), "remapped", len(other), "not found", len(bad), "multiple names", len(dupe)
print
print "DUPE:", dupe
print
print "BAD:", bad
print
print "OTHER:", other

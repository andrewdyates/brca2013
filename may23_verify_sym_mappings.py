#!/usr/bin/python
"""
good: 13632 remapped 187 not found 396 multiple names 3
"""
import hugo_gene_symbols
H = hugo_gene_symbols.load()
FNAME = "E7307.original.row.names.txt"

good = set()
other = {}
bad = set()
dupe = {}
for line in open(FNAME):
  s = line.strip()
  z = H.find_sym(s, allow_dupe=True)
  if s != z:
    if isinstance(z,set):
      dupe[s] = z
    elif z:
      other[s] = z
    else:
      bad.add(s)
  else:
    good.add(s)

print "good:", len(good), "remapped", len(other), "not found", len(bad), "multiple names", len(dupe)
print
print "DUPE:", dupe
print
print "BAD:", bad
print
print "OTHER:", other

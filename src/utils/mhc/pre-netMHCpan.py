#!/usr/bin/env python3
import os, sys

db = sys.argv[1]
hlastr = sys.argv[2]

dbdict = {}
dbfile = open(db, 'r')
if dbfile:
    for line in dbfile.readlines():
        arr = line.strip().split()
        dbdict[arr[0]] = arr[1]
    dbfile.close()
else:
    sys.stderr.write("Cannot open %s\n" % db)
    sys.exit(1)

hlist = hlastr.strip().split(",")
out = []
for hla in hlist:
    if hla in dbdict:
        out.append(hla)
print(",".join(out))

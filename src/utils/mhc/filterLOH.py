#!/usr/bin/env python3
import sys

HLALOH = sys.argv[1]
neofile = sys.argv[2]

loh = []
fd = open(HLALOH, "r", encoding='utf8')
for line in fd.readlines():
    arr = line.strip().split()
    if arr[1] != "阳性":
        continue
    hla = arr[0].split("_")
    hlastr = hla[0].upper() + "-" + hla[1].upper() + "*" + hla[2] + ":" + hla[3]
    loh.append(hlastr)
fd.close()

fd = open(neofile, "r", encoding='utf8')
for line in fd.readlines():
    arr = line.strip().split("\t")
    if arr[2] in loh:
        continue
    print(line.strip())
fd.close()

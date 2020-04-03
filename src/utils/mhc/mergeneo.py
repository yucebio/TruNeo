#!/usr/bin/env python3
import sys

mutneof = sys.argv[1]
fusneof = sys.argv[2]

neoarr = []
fd = open(fusneof, "r")
for line in fd.readlines():
    arr = line.strip().split("\t")
    if arr[0] == "chrom1":
        continue
    record = {}
    record["HLA"] = arr[15]
    record["Peptide"] = arr[13]
    record["NeoepitopeScore"] = arr[14]
    if arr[14] != "NA" and float(arr[14]) <= 50:
        record["NeoepitopeType"] = "SB"
    else:
        record["NeoepitopeType"] = "WB"
    record["Gene"] = arr[6]
    record["Psite"] = arr[16]
    neoarr.append(record)
fd.close()

fd = open(mutneof, "r")
head = []
glist = []
lines = fd.readlines()
if len(lines) > 1:
    count = int(lines[-1].strip().split("\t")[0])
else:
    count = 0
for line in lines:
    arr = line.strip().split("\t")
    if arr[0] == "NeoRank":
        head = arr
    print(line.strip())
    if arr[8] not in glist:
        glist.append(arr[8])
fd.close()


for i in neoarr:
    if i["Gene"] in glist:
        continue
    count += 1
    outmp = []
    for j in head:
        if j == "NeoRank":
            outmp.append(str(count))
        elif j in i:
            outmp.append(i[j])
        else:
            outmp.append("NA")
    print("\t".join(outmp))

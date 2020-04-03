#!/usr/bin/env python3
import sys
import re

sortneo = sys.argv[1]
detailf = sys.argv[2]

peptide = {}
mutinfo = {}
fd = open(sortneo, "r", encoding='utf8')
num = 0
for line in fd.readlines():
    arr = line.strip().split("\t")
    if arr[0] == "NeoRank":
        continue
    if float(arr[2]) == 0:
        continue
    num += 1
    keys = arr[12] + "\t" + arr[13] + "\t" + arr[14] + "\t" + arr[15] + "\t" + arr[19].replace("p.", "") + "\t" + arr[5] + "\t" + arr[6]
    mutinfo[arr[0]] = {}
    if not arr[7] in peptide:
        peptide[arr[7]] = []
    peptide[arr[7]].append(arr[0])
    mutinfo[arr[0]]["s"] = (int(arr[4]) - 1) * 3
    mutinfo[arr[0]]["e"] = (int(arr[4]) + len(arr[6]) - 1) * 3
    mutinfo[arr[0]]["i"] = keys
fd.close()

fd = open(detailf, "r", encoding='utf8')
pID = "na"
seq = ""
readf = 0
for line in fd.readlines():
    if re.match(">", line.strip()):
        arr = line.strip().split()
        if arr[4] == "MTsubcdna":
            readf = 1
        else:
            readf = 0
            continue

        if pID in peptide:
            for i in peptide[pID]:
                mutinfo[i]["i"] += "\t" + seq[mutinfo[i]["s"]:mutinfo[i]["e"]]

        pID = arr[0].replace(">", "")
        seq = ""
    elif readf == 1:
        seq += line.strip()
fd.close()
if pID in peptide:
    for i in peptide[pID]:
        mutinfo[i]["i"] += "\t" + seq[mutinfo[i]["s"]:mutinfo[i]["e"]]

print("Gene\tChr\tPosition\tTranscript\tProteinChangeAbbr\tHLA\tPeptide\tcdna")
for i in sorted(mutinfo.keys(), key = lambda x:int(x)):
    print(mutinfo[i]["i"])

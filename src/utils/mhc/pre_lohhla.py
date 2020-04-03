#!/usr/bin/env python3
import sys

sample = sys.argv[1]
hlafile = sys.argv[2]
ascat = sys.argv[3]
outd = sys.argv[4]

HLA = []
with open(hlafile, "r") as fd:
    for line in fd.readlines():
        arr = line.strip().split()
        HLA.append(arr[1])
        HLA.append(arr[2])
set1 = set(HLA)
out1 = open(outd + "/hla.list", "w")
sort1 = sorted(set1)
out1.write("{}\n".format("\n".join(sort1)))

out2 = open(outd + "/CN.solution", "w")
out2.write("Ploidy\ttumorPurity\ttumorPloidy\n")
with open(ascat, "r") as fd:
    for line in fd.readlines():
        arr = line.strip().split()
        if arr[0] == "Ploidy":
            ploidy = arr[1]
        if arr[0] == "rho":
            rho = arr[1]
if rho == "1":
    rho = "0.2"
out2.write("{}\t{}\t{}\n".format(sample, rho, ploidy))

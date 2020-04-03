#!/usr/bin/env python3
import os, sys

infile = sys.argv[1]
with open(infile, "r") as fd:
    event = []
    for line in fd.readlines():
        arr = line.strip().split("\t")
        if arr[0] == "#FusionName":
            sys.stderr.write(line)
            continue
        key = arr[0].split("--")[1] + "--" + arr[0].split("--")[0]
        if key in event or arr[0] in event:
            continue
        if arr[3] == "ONLY_REF_SPLICE" and int(arr[1]) < 5:
            continue
        if arr[3] == "INCL_NON_REF_SPLICE" and (int(arr[1]) < 5 or int(arr[2]) < 2 or int(arr[1]) + int(arr[2]) < 10):
            continue
        info1 = arr[5].split(":")
        info2 = arr[7].split(":")
        if info1[0] == info2[0] and abs(int(info1[1]) - int(info2[1])) < 10000000:
            continue
        event.append(arr[0])
        sys.stderr.write(line)
        name = arr[0].replace("--", ">>")
        chr1 = info1[0].replace("chr", "")
        chr2 = info2[0].replace("chr", "")
        if info1[2] == "+" and info2[2] == "+":
            print("\t".join((chr1, "-1", info1[1], chr2, str(int(info2[1]) - 1), "-1", name, str(int(arr[1]) + int(arr[2])), "+", "+", arr[9])))
        elif info1[2] == "-" and info2[2] == "-":
            print("\t".join((chr1, str(int(info1[1]) - 1), "-1", chr2, "-1", info2[1], name, str(int(arr[1]) + int(arr[2])), "-", "-", arr[9])))
        elif info1[2] == "+" and info2[2] == "-":
            print("\t".join((chr1, "-1", info1[1], chr2, "-1", info2[1], name, str(int(arr[1]) + int(arr[2])), "+", "-", arr[9])))
        else:
            print("\t".join((chr1, str(int(info1[1]) - 1), "-1", chr2, str(int(info2[1]) - 1), "-1", name, str(int(arr[1]) + int(arr[2])), "-", "+", arr[9])))

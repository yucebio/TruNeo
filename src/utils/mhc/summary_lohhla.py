#!/usr/bin/env python3
import sys, os

hlaloh = sys.argv[1]
hlalist = sys.argv[2]

positive = []
if os.path.exists(hlaloh):
    with open(hlaloh, "r") as fd:
        for line in fd.readlines():
            arr = line.strip().split()
            if arr[0] == "region" or arr[-6] == "NA" or float(arr[-6]) >= 0.0002:
                continue
            positive.append(arr[-4])

print("血液样本HLA分型\t肿瘤样本HLA缺失检测结果")
with open(hlalist, "r") as fd:
    for line in fd.readlines():
        htype = line.strip()
        if htype in positive:
            judge = "阳性"
        else:
            judge = "阴性"
        print("{}\t{}".format(htype, judge))

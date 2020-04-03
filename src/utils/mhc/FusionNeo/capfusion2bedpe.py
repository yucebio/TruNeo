#!/usr/bin/env python3
import os, sys

def annoexon(ref, gene, pos, orien):
    posout = {}
    with open(ref, "r") as fd:
        for line in fd.readlines():
            arr = line.strip().split("\t")
            if arr[12] != gene:
                continue
            if int(arr[4]) >= int(pos) or int(arr[5]) < int(pos):
                continue
            starts = arr[9].strip(",").split(",")
            ends = arr[10].strip(",").split(",")
            for i in range(len(starts)):
                if int(pos) > int(starts[i]) and int(pos) <= int(ends[i]):
                    if orien == "l":
                        posout[pos] = 1
                    else:
                        posout[str(int(pos)-1)] = 1
                    break
                if int(pos) < int(ends[i]):
                    if orien == "l":
                        posout[ends[i-1]] = 1
                    else:
                        posout[starts[i]] = 1
                    break
    if len(posout) == 0:
        if orien == "l":
            posout[pos] = 1
        else:
            posout[str(int(pos)-1)] = 1
    return posout

def main():
    infile = sys.argv[1]
    refgene = sys.argv[2]
    with open(infile, "r") as fd:
        for line in fd.readlines():
            arr = line.strip().split("\t")
            if arr[0] == "Gene1":
                if arr[10] == "Freq(%)":
                    frei = 10
                else:
                    frei = 8
                continue
            if arr[-1] == "Unknown":
                continue
            genes = arr[-1].split(">>")
            if genes[0] == arr[0]:
                chr1 = arr[1].replace("chr", "")
                pos1 = arr[2]
                strand1 = arr[-2][0]
                chr2 = arr[4].replace("chr", "")
                pos2 = arr[5]
                strand2 = arr[-2][1]
            else:
                chr1 = arr[4].replace("chr", "")
                pos1 = arr[5]
                strand1 = arr[-2][1]
                chr2 = arr[1].replace("chr", "")
                pos2 = arr[2]
                strand2 = arr[-2][0]
            if strand1 == "+" and strand2 == "+":
                posarr1 = annoexon(refgene, genes[0], pos1, "l")
                posarr2 = annoexon(refgene, genes[1], pos2, "r")
                for p1 in posarr1:
                    for p2 in posarr2:
                        print("\t".join((chr1, "-1", p1, chr2, p2, "-1", arr[-1], arr[6], strand1, strand2, arr[frei])))
            elif strand1 == "-" and strand2 == "-":
                posarr1 = annoexon(refgene, genes[0], pos1, "r")
                posarr2 = annoexon(refgene, genes[1], pos2, "l")
                for p1 in posarr1:
                    for p2 in posarr2:
                        print("\t".join((chr1, p1, "-1", chr2, "-1", p2, arr[-1], arr[6], strand1, strand2, arr[frei])))
            elif strand1 == "+" and strand2 == "-":
                posarr1 = annoexon(refgene, genes[0], pos1, "l")
                posarr2 = annoexon(refgene, genes[1], pos2, "l")
                for p1 in posarr1:
                    for p2 in posarr2:
                        print("\t".join((chr1, "-1", p1, chr2, "-1", p2, arr[-1], arr[6], strand1, strand2, arr[frei])))
            else:
                posarr1 = annoexon(refgene, genes[0], pos1, "r")
                posarr2 = annoexon(refgene, genes[1], pos2, "r")
                for p1 in posarr1:
                    for p2 in posarr2:
                        print("\t".join((chr1, p1, "-1", chr2, p2, "-1", arr[-1], arr[6], strand1, strand2, arr[frei])))

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import sys
import re
import os

def usage():
    print("Convert transvar output to netmhcpan input.")
    print("Usage: convert.py <input_file> <output_prefix> [output_path] [debug]")

def getid(linelist, pepno, compre3):
    outline = [str(pepno), linelist[2], '', '', '', '', '_', '_', '', '_']
    m = re.search(r'^(\S+)', linelist[1])
    if m:
        outline[4] = m.group(1)
    m = compre3.search(linelist[0])
    if m:
        outline[2] = m.group(1)
        outline[3] = m.group(2)
    poslist = linelist[4].split('/')
    if len(poslist) == 3:
        outline[5] = poslist[1]
        outline[8] = poslist[2]
    return "\t".join(outline)


if len(sys.argv) < 2:
    usage()
    sys.exit(1)
elif len(sys.argv) < 3:
    path = '.'
else:
    path = sys.argv[3]

infile = sys.argv[1]
prefix = sys.argv[2]
wildfile = path + '/' + prefix + ".widetype.txt"
wildfd = open(wildfile, "w")
mutfile = path + '/' + prefix + ".result.txt"
mutfd = open(mutfile, 'w')
idfile = path + '/' + "peptide.peptideID.txt"
idfd = open(idfile, 'w')
failfile = path + '/' + prefix + ".failed.txt"
failfd = open(failfile, 'w')

pepno = 1
compre1 = re.compile(r'(\w+)__\[frameshift_([^\>]+)\>([^\]]+)\]$')
compre2 = re.compile(r'(\w*)__\[([^\>]+)\>([^\]]+)\]_*(\S*)$')
compre3 = re.compile(r'chr(\w+):g\.(\d+)')
compre4 = re.compile(r'p\.(\w)\d+\*')
with open(infile, 'r') as fd:
    for line in fd:
        linelist = line.strip().split("\t")
        if linelist[0] == "input":
            continue
        m = re.search(r'variant_protein_seq=([^;]+);', linelist[6])
        if not m:
            failfd.write(line)
            continue
        mm = m.group(1)
        if mm.find("frameshift") > -1:
            m = compre1.search(mm)
        else:
            m = compre2.search(mm)
        if m:
            leftstr = m.group(1)
            wildstr = m.group(2)
            mutstr = m.group(3)
            rightstr = ''
            if len(m.groups()) > 3:
                rightstr = m.group(4)
            if len(leftstr) >= 13:
                leftstr = leftstr[-13:]
            if len(rightstr) >= 13:
                rightstr = rightstr[:13]
            mm = compre4.search(linelist[4])
            if mm:
                wildstr = mm.group(1) + wildstr
            if len(wildstr) >= 13:
                wildstr = wildstr[:13]
                wildpep = leftstr + wildstr
            else:
                wildpep = leftstr + wildstr + rightstr
            if len(mutstr) >= 13:
                mutstr = mutstr[:13]
                mutpep = leftstr + mutstr
            else:
                mutpep = leftstr + mutstr + rightstr
            wildpep = wildpep.replace('*', '')
            mutpep = mutpep.replace('*', '')
            
            pepid = ">Peptide" + str(pepno)
            wildfd.write(pepid + "\n")
            wildfd.write(wildpep + "\n")
            mutfd.write(pepid + "\n")
            mutfd.write(mutpep + "\n")
            idline = getid(linelist, pepno, compre3)
            idfd.write(idline + "\n")
            if len(sys.argv) > 3:
                print(idline)
                print(wildpep)
                print(mutpep)

            pepno = pepno + 1
        else:
            failfd.write(mm + "\n")

wildfd.close()
mutfd.close()
idfd.close()
failfd.close()

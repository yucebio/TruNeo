#!/usr/bin/env python3
import os, sys
import getopt
import pysam
import re

def usage():
    print("Do the snv phasing")
    print("Author: chenly")
    print("Usage: snv_phasing.py [options]")
    print("")
    print(" -i <file>    Input snv result file, sorted by coordinate. Required")
    print(" -b <file>    Input sample bam file. Required")
    print(" -a <file>    Input snv annodb file. Required")
    print(" -x <file>    Input indel result file. Optional")
    print(" -c <int>     Index num of Protein in annodb file. Default 12")
    print(" -r <file>    Input reference fasta file. Required")
    print(" -d <int>     Max distance for phasing. Default 45")
    print(" -t <str>     Input file type. Could be varscan. Default varscan")
    print(" -f <float>   Min variant allele frequency. Default 0.03")
    print(" -o <str>     Output prefix")

def get_qname(bamh, reg, pos, alt, end):
    altlist = []
    otherlist = []
    for column in bamh.pileup(region=reg):
        if column.pos == pos:
            for reads in column.pileups:
                qpos = reads.query_position
                align = reads.alignment
                if bin(align.flag)[-1] == "1" and not reads.is_del and not reads.is_refskip and reads.indel == 0 and not reads.is_head and not reads.is_tail and align.query_qualities[qpos] >= 20 and align.query_length >= 70 and not re.search("H", align.cigarstring) and len(re.findall("I|D", align.cigarstring)) < 3 and align.reference_name == align.next_reference_name and abs(align.reference_start - align.next_reference_start) < 1000 and ((int(end) > pos and align.reference_end >= int(end)) or (int(end) <= pos and align.reference_start < int(end))):
                    if align.query_sequence[qpos] == alt:
                        altlist.append(align.query_name)
                    else:
                        otherlist.append(align.query_name)
    set1 = set(altlist)
    set2 = set(otherlist)
    return set1 - set2

def get_qname2(bamh, reg, pos, alt, end):
    altlist = []
    for column in bamh.pileup(region=reg):
        if column.pos == pos:
            for reads in column.pileups:
                qpos = reads.query_position
                align = reads.alignment
                if bin(align.flag)[-1] == "1" and not reads.is_refskip and reads.indel == alt and not reads.is_head and not reads.is_tail and align.query_qualities[qpos] >= 20 and align.query_length >= 70 and not re.search("H", align.cigarstring) and len(re.findall("I|D", align.cigarstring)) < 3 and align.reference_name == align.next_reference_name and abs(align.reference_start - align.next_reference_start) < 1000 and ((int(end) > pos and align.reference_end >= int(end)) or (int(end) <= pos and align.reference_start < int(end))):
                    altlist.append(align.query_name)
    set1 = set(altlist)
    return set1

def phasing(set1, set2):
    overlap = len(set1 & set2)
    if len(set1) > 0 and overlap / len(set1) >= 0.75:
        return 2, overlap
    elif len(set2) > 0 and overlap / len(set2) >= 0.75:
        return 1, overlap
    else:
        return 0, overlap

def read_anno(inf, index):
    dictmp = {}
    for line in inf.readlines():
        arr = line.strip().split("\t")
        if arr[0] == 'Chr':
            continue
        arr[0] = 'chr' + arr[0]
        key = arr[0] + ':' + arr[1]
        m = re.search("p\.(\D+\d+)\D+", arr[int(index)])
        if m:
            dictmp[key] = m.group(1)
        else:
            dictmp[key] = 'na'
    return dictmp

def readvar(inf, idict):
    header = ""
    with open(inf, "r") as fd:
        for line in fd.readlines():
            arr = line.strip().split("\t")
            if arr[0] == "chrom":
                header = line
                continue
            if re.match("\+", arr[3]):
                arr[3] = arr[2] + arr[3].replace("+", "")
                arr[11] = arr[2] + "/" + arr[3]
            elif re.match("\-", arr[3]):
                arr[2] = arr[2] + arr[3].replace("-", "")
                arr[3] = arr[2][0]
                arr[11] = arr[2] + "/" + arr[3]
            linfo = "\t".join(arr)
            if not arr[0] in idict:
                idict[arr[0]] = {}
            idict[arr[0]][arr[1]] = linfo + "\n"
    return header

def read_varscan(inf, bamh, fah, dist, minfre, outp, annodict, idict):
    chrom = 'na'
    postmp = 0
    infotmp = 'na'
    outtmp = 'na'
    cnum = 0
    cflag = 0
    cchr = 'na'
    cpos = ''
    calt = ''
    out1 = open(outp + '.snv.new.xls', 'w')
    out2 = open(outp + '.phasing.info', 'w')
    if inf == outp + '.merge.temp':
        out3 = open(outp + '.indel.new.xls', 'w')
    infile = open(inf, "r")
    for line in infile.readlines():
        arr = line.strip().split("\t")
        if arr[0] == 'chrom':
            out1.write(line)
            if inf == outp + '.merge.temp':
                out3.write(line)
            continue
        reflen = 0
        arrtmp = []
        if infotmp != 'na':
            arrtmp = infotmp.split("\t")
            if len(arrtmp[3]) < len(arrtmp[2]):
                reflen = len(arrtmp[2]) - len(arrtmp[3])
        diff = int(arr[1]) - int(postmp) - reflen + 1
        if chrom == 'na':
            outtmp = line
        elif chrom != arr[0] or (chrom == arr[0] and (diff > dist or diff <= 0)):
            arrout = outtmp.strip().split("\t")
            if len(arrout[2]) == len(arrout[3]):
                out1.write(outtmp)
            else:
                out3.write(outtmp)
            if cflag == 1:
                out2.write("%d\t%s\t%s\t%s\n" % (cnum, cchr, cpos, calt))
                cflag = 0
            outtmp = line
        elif chrom == arr[0] and diff <= dist:
            key1 = chrom + ':' + postmp
            key2 = arr[0] + ':' + arr[1]
            region = chrom + ':' + postmp + '-' + postmp
            if len(arrtmp[2]) == len(arrtmp[3]):
                set1 = get_qname(bamh, region, int(postmp) - 1, arrtmp[3], arr[1])
            else:
                indlen = len(arrtmp[3]) - len(arrtmp[2])
                set1 = get_qname2(bamh, region, int(postmp) - 1, indlen, arr[1])
            region = arr[0] + ':' + arr[1] + '-' + arr[1]
            if len(arr[2]) == len(arr[3]):
                set2 = get_qname(bamh, region, int(arr[1]) - 1, arr[3], postmp)
            else:
                indlen = len(arr[3]) - len(arr[2])
                set2 = get_qname2(bamh, region, int(arr[1]) - 1, indlen, postmp)
            flag, overlap = phasing(set1, set2)
            outarr = outtmp.strip().split("\t")
            if (int(arr[1]) - int(postmp) - reflen < 9 and (len(outarr[2]) != len(outarr[3]) or len(arr[2]) != len(arr[3]))) or (len(outarr[2]) == len(outarr[3]) and len(arr[2]) == len(arr[3]) and key1 in annodict and key2 in annodict and annodict[key1] == annodict[key2] and annodict[key1] != 'na'):
                if flag != 0:
                    if int(arr[1]) - 1 > int(postmp) + reflen:
                        ref = fah.fetch(arr[0], int(postmp) + reflen, int(arr[1]) - 1)
                        outarr[2] = outarr[2] + ref.upper() + arr[2]
                        outarr[3] = outarr[3] + ref.upper() + arr[3]
                        outarr[7] = outarr[7] + ref.upper() + arr[7]
                    elif len(outarr[2]) == len(outarr[3]) and len(arr[2]) == len(arr[3]):
                        outarr[2] = outarr[2] + arr[2]
                        outarr[3] = outarr[3] + arr[3]
                        outarr[7] = outarr[7] + arr[7]
                    if flag == 1:
                        outarr[4] = arr[4]
                        outarr[5] = arr[5]
                        outarr[6] = arr[6]
                        outarr[14] = arr[14]
                    outarr[8] = (int(outarr[8]) + int(arr[8])) // 2
                    outarr[9] = overlap
                    outarr[10] = str('{:.2f}'.format(100 * outarr[9] / (outarr[8] + outarr[9]))) + '%'
                    if len(outarr[2]) == len(outarr[3]) and len(arr[2]) == len(arr[3]):
                        if int(arr[1]) - 1 > int(postmp) + reflen:
                            outarr[11] = outarr[11] + ref.upper() + arr[11]
                        else:
                            outarr[11] = outarr[11] + arr[11]
                    else:
                        outarr[11] = outarr[2] + "/" + outarr[3]
                    for i, v in enumerate(outarr):
                        outarr[i] = str(v)
                    outtmp = '\t'.join(outarr) + '\n'

                    if cflag == 0:
                        cnum += 1
                        cchr = chrom
                        cpos = outarr[1]
                        calt = outarr[3]
                        cflag = 1
                    elif int(arr[1]) - 1 <= int(postmp) + reflen:
                        calt = calt + arr[3]
                    else:
                        calt = calt + ref.upper() + arr[3]
                else:
                    if len(outarr[2]) == len(outarr[3]):
                        out1.write(outtmp)
                    else:
                        out3.write(outtmp)
                    if cflag == 1:
                        out2.write("%d\t%s\t%s\t%s\n" % (cnum, cchr, cpos, calt))
                        cflag = 0
                    outtmp = line

                if flag == 1:
                    remain = len(set1 - set2)
                    ff = rest_judge(remain, arrtmp, minfre)
                    if ff != '0' and overlap / len(set2) <= 0.9:
                        if len(arrtmp[2]) == len(arrtmp[3]):
                            out1.write(ff)
                        else:
                            out3.write(ff)
                elif flag == 2:
                    remain = len(set2 - set1)
                    ff = rest_judge(remain, arr, minfre)
                    if ff != '0' and overlap / len(set1) <= 0.9:
                        if len(arr[2]) == len(arr[3]):
                            out1.write(ff)
                        else:
                            out3.write(ff)
            elif len(outarr[2]) != len(outarr[3]) or len(arr[2]) != len(arr[3]) or (len(outarr[2]) == len(outarr[3]) and len(arr[2]) == len(arr[3]) and ((key1 in annodict and annodict[key1] != 'na') or (key2 in annodict and annodict[key2] != 'na'))):
                if flag != 0:
                    if cflag == 0:
                        cnum += 1
                        cchr = chrom
                        cpos = postmp + ',' + arr[1]
                        calt = arrtmp[3] + ',' + arr[3]
                        cflag = 1
                    else:
                        cpos = cpos + ',' + arr[1]
                        calt = calt + ',' + arr[3]
                elif cflag == 1:
                    out2.write("%d\t%s\t%s\t%s\n" % (cnum, cchr, cpos, calt))
                    cflag = 0
                if len(outarr[2]) == len(outarr[3]):
                    out1.write(outtmp)
                else:
                    out3.write(outtmp)
                outtmp = line
            else:
                out1.write(outtmp)
                if cflag == 1:
                    out2.write("%d\t%s\t%s\t%s\n" % (cnum, cchr, cpos, calt))
                    cflag = 0
                outtmp = line
        chrom = arr[0]
        postmp = arr[1]
        infotmp = line.strip()
    if outtmp != 'na':
        arrout = outtmp.strip().split("\t")
        if len(arrout[2]) == len(arrout[3]):
            out1.write(outtmp)
        else:
            out3.write(outtmp)
    if cflag == 1:
        out2.write("%d\t%s\t%s\t%s\n" % (cnum, cchr, cpos, calt))
        cflag = 0
    infile.close()
    out1.close()
    out2.close()
    if inf == outp + '.merge.temp':
        out3.close()

def rest_judge(altn, array, minfre):
    if int(array[8]) + altn > 9 and altn > 2:
        fre = altn / (int(array[8]) + altn)
        if fre >= minfre:
            array[9] = altn
            array[10] = str('{:.2f}'.format(100 * fre)) + '%'
            for i, v in enumerate(array):
                array[i] = str(v)
            return '\t'.join(array) + '\n'
        else:
            return '0'
    else:
        return '0'

def main():
    if len(sys.argv) < 9:
        usage()
        sys.exit(1)

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:b:a:x:c:r:d:t:f:o:')
    except getopt.GetoptError:
        sys.exit(1)

    config = {}
    config['pindex'] = 12
    config['dist'] = 45
    config['type'] = 'varscan'
    config['vaf'] = 0.03
    for opt, value in opts:
        if opt == '-i':
            config['infile'] = value
        elif opt == '-b':
            config['bam'] = value
        elif opt == '-a':
            config['anno'] = value
        elif opt == '-x':
            config['indel'] = value
        elif opt == '-c':
            if re.search("\D", value):
                usage()
                sys.exit(1)
            config['pindex'] = value
        elif opt == '-r':
            config['ref'] = value
        elif opt == '-d':
            if re.search("\D", value):
                usage()
                sys.exit(1)
            config['dist'] = int(value)
        elif opt == '-t':
            config['type'] = value
        elif opt == '-f':
            if re.search("[A-Za-z_]", value):
                usage()
                sys.exit(1)
            config['vaf'] = float(value)
        elif opt == '-o':
            config['outp'] = os.path.abspath(value)
    if 'bam' not in config or 'anno' not in config or 'ref' not in config:
        usage()
        sys.exit(1)

    annofile = open(config['anno'], 'r', encoding='utf8')
    annodb = {}
    if annofile:
        annodb = read_anno(annofile, config['pindex'])
        annofile.close()
    else:
        sys.stderr.write("Cannot open %s\n" % config['anno'])
        sys.exit(1)

    snvfile = config['infile']
    idict = {}
    if 'indel' in config:
        head = readvar(config['infile'], idict)
        readvar(config['indel'], idict)
        outmerge = open(config['outp'] + '.merge.temp', 'w')
        outmerge.write(head)
        for i in sorted(idict.keys()):
            for j in sorted(idict[i].keys(), key = lambda x:int(x)):
                outmerge.write(idict[i][j])
        outmerge.close()
        snvfile = config['outp'] + '.merge.temp'

    bamfile = pysam.AlignmentFile(config['bam'], 'rb')
    fafile = pysam.FastaFile(config['ref'])
    read_varscan(snvfile, bamfile, fafile, config['dist'], config['vaf'], config['outp'], annodb, idict)
    bamfile.close()
    fafile.close()
    if 'indel' in config:
        os.remove(snvfile)

if __name__ == '__main__':
    main()

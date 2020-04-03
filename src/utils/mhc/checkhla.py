#!/usr/bin/env python3

import sys
#import json
import os
from collections import defaultdict
import getopt

"""
如果 polysolver 预测的 tumor 和 normal 一致，输出 polysolver 结果。
如果 polysolver 预测的 tumor 是 normal 的子集，输出 polysolver 的 normal 结果。
如果 polysolver 预测的 tumor 和 normal 不一致，tumor 也不是 normal 的子集，就看 bwahla 预测结果，如果 bwahla 预测的 tumor 和 normal 一致就输出 bwahla 结果，如果 bwahla 预测的 tumor 和 normal 不一致，就不输出预测结果。
"""
def usage():
	print("Usage: /usr/bin/python3.4 %s [options]"%(sys.argv[0]))
	print("")
	print("Options:")
	print("")
	print("  -d    hla typing result folder, this folder contains polysolver and bwahla outputdir")
	print("  -t    tumor sample name")
	print("  -b    control sample name")
	print("  -s    somatic symbol")
	print("  -o    output dir, this script will generate a text file like: ")
	print("        outputdir/somatic.hla1.result")
	print("  -h    print usage")

if len(sys.argv) == 1:
	usage()
	sys.exit(1)

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:t:b:s:o:h")
except getopt.GetoptError as err:
	print(err)
	usage()
	sys.exit(1)

options = {}
for option,value in opts:
	if option == "-d":
		options["hladir"] = value
	elif option == "-t":
		options["tumor"] = value
	elif option == "-b":
		options["control"] = value
	elif option == "-s":
		options["somatic"] = value
	elif option == "-o":
		options["outputdir"] = value
	elif option == "-h":
		usage()
		sys.exit()
	else:
		assert False, "unhandled option"



"""
with open(json_file, 'r') as proj_config_file:
	proj_config = json.load(proj_config_file)
	samplelist = []
	workdir = proj_config["project"]["workdir"]
	for samplepair in proj_config["somatic"]:
		samplelist.append(samplepair.split("-VS-"))
"""

def polysolverparser(hlafile):
	polysolver_hlatype = defaultdict(set)
	with open(hlafile, 'r') as hla_result:
		for line in hla_result:
			hla_class, *hla_type = line.strip().split()
			for hla_subtype in hla_type:
				hla_subtype = hla_subtype.lstrip("hlabc_")
				hla_subtype = ":".join(hla_subtype.split("_")[:2])
				polysolver_hlatype[hla_class].add(hla_subtype)
	
	return polysolver_hlatype

def bwahlaparser(hlafile):
	bwahla_hlatype = defaultdict(set)
	with open(hlafile, 'r') as hla_result:
		for line in hla_result:
			prefix, hla_subtype1, hla_subtype2, *stats = line.strip().split()
			hla_class = hla_subtype1.split("*")[0]
			hla_subtype1 = hla_subtype1.split("*")[1].split(":")[:2]
			hla_subtype1 = ":".join(hla_subtype1)
			hla_subtype2 = hla_subtype2.split("*")[1].split(":")[:2]
			hla_subtype2 = ":".join(hla_subtype2)

			bwahla_hlatype[hla_class].add(hla_subtype1)
			bwahla_hlatype[hla_class].add(hla_subtype2)

	return bwahla_hlatype

"""
for pair in samplelist:
	tumor = pair[0]
	control = pair[1]
"""	
hladir = options["hladir"]
tumor = options["tumor"]
control = options["control"]
outputdir = options["outputdir"]

polysolver_dir = hladir + "/polysolver/"
bwahla_dir     = hladir + "/bwahla/"

tumor_hla_polysolver = polysolverparser(polysolver_dir + "/" + tumor + "/" + tumor + ".hla.txt")
control_hla_polysolver = polysolverparser(polysolver_dir + "/" + control + "/" + control + ".hla.txt")

tumor_hla_bwahla = bwahlaparser(bwahla_dir + "/" + tumor + "/" + tumor + ".hla.top")
control_hla_bwahla = bwahlaparser(bwahla_dir + "/" + control + "/" + control + ".hla.top")

pair_result_prefix = options["somatic"]
pair_result_dir = outputdir

final_hlatype_filename = pair_result_dir + "/" + pair_result_prefix + ".hla1.result"
# check path exists or not
if os.path.exists(pair_result_dir):
	final_hlatype = open(final_hlatype_filename, 'w')
else:
	os.makedirs(pair_result_dir)
	final_hlatype = open(final_hlatype_filename, 'w')
#final_hlatype_filename = pair_result_dir + "/" + pair_result_prefix + ".hla1.result"
#final_hlatype = open(final_hlatype_filename, 'w')

final_hlatype_content = []
for hla_class in ["HLA-A","HLA-B","HLA-C"]:
	if tumor_hla_polysolver[hla_class] == control_hla_polysolver[hla_class]:
		final_hlatype_content.extend([hla_class + i for i in  list(tumor_hla_polysolver[hla_class])])

	elif tumor_hla_polysolver[hla_class] < control_hla_polysolver[hla_class]:
                final_hlatype_content.extend([hla_class + i for i in  list(control_hla_polysolver[hla_class])])
                '''
		diff_hla_polysolver = control_hla_polysolver[hla_class] - tumor_hla_polysolver[hla_class]
		if diff_hla_polysolver < control_hla_bwahla[hla_class] or diff_hla_polysolver < tumor_hla_bwahla[hla_class]:
			final_hlatype_content.extend([hla_class + i for i in  list(control_hla_polysolver[hla_class])])
		else:
			if tumor_hla_bwahla[hla_class] == control_hla_bwahla[hla_class]:
				final_hlatype_content.extend([hla_class + i for i in  list(control_hla_bwahla[hla_class])])
                '''
	else:
		if tumor_hla_bwahla[hla_class] == control_hla_bwahla[hla_class]:
			final_hlatype_content.extend([hla_class + i for i in  list(control_hla_bwahla[hla_class])])

final_hlatype.write(",".join(final_hlatype_content) + "\n")
final_hlatype.close()

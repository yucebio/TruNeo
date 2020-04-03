#=============================================
#
# File: run.TruNeo.single.sh
#
# Version: 1.0
#
# Original Author: Jiaqian Wang, wangjq@yucebio.com
#
# Organization: Yucebio Inc., Shenzhen
#
# Description: This is the main scripts of TruNeo, 
# which aimed at predicting and prioritizing 
# neoantigens from DNA Variation & RNA-Seq.
#
#=============================================

#1.Prerequisite
#1) Please install the following softwares, and configure them properly,
#	You may need to configure the path in src/config/*:
#	·Perl 5.22.1
#	·R 3.4.4
#	·python=Python3.6
#	·python2=Python2.7
#	·netMHCpan3.2
#		set PROG_NETMHC=/path/to/netMHCpan
#	·netMHC 4.0
#	·netchop 3.1
#		set UTIL_TAP_PREDICT=/path/to/predict.py
#	·Star-fusion
#	·HLAminer v1.3
#	
#2) The freely available resources, and could be helpful to 
#   generate inputs. You may need to configure the path in src/config/*:
#	·Reference genome hg19, download from UCSC
#	·human_g1k_v37_decoy.fasta, downloaded from GATK
#	·Annotation Databases, donwload from GATK and GENCODE
#	·STAR and Star-fusion official reference
#
#3) The format of in-house databases are provided as-is in ./src/yucemed, 
#	or set the following paths in src/config/*:
#
#	DB_NONCHR_FA=/path/to/human_g1k_v37_decoy.fasta
#	DB_TRANS_FASTA=/path/to/hg19.gencodev27.transcripts.fa
#	DB_GENEPRED=/path/to/hg19.gencodev27.protein_coding.genePred
#	DB_REFGENE=/path/to/refGene.hg19.gencodev27.txt
#
#2.Usage
#	bash run.TruNeo.single.sh
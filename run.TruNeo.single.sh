#!/bin/bash
#=============================================
# File: run.TruNeo.single.sh
# Version: 1.0
# Original Author: Jiaqian Wang, wangjq@yucebio.com
# Organization: Yucebio Inc., Shenzhen
# Description: This is the main scripts of TruNeo, which aimed at predicting and prioritizing neoantigens from DNA Variation & RNA-Seq.
#=============================================

#1.Prerequisite
#1) Please install the following softwares, and configure them properly:
#	Perl 5.22.1
#	R 3.4.4
#	python=Python3.6
#	python2=Python2.7
#	netMHCpan3.2
#		set PROG_NETMHC=/path/to/netMHCpan
#	netMHC 4.0
#	netchop 3.1
#		set UTIL_TAP_PREDICT=/path/to/predict.py
#	Star-fusion
#	HLAminer v1.3
#	
#2)	The freely available resources, and could be helpful to generate inputs. You may need to configure the path in src/config/*:
#	Reference genome hg19, download from UCSC
#	human_g1k_v37_decoy.fasta, downloaded from GATK
#	Annotation Databases, donwload from GATK and GENCODE
#	STAR and Star-fusion official reference

#3) Set paths to the following databases:
#	DB_NONCHR_FA=/path/to/human_g1k_v37_decoy.fasta
#	DB_TRANS_FASTA=/path/to/hg19.gencodev27.transcripts.fa
#	DB_GENEPRED=/path/to/hg19.gencodev27.protein_coding.genePred
#	DB_REFGENE=/path/to/refGene.hg19.gencodev27.txt

#2.ENV
#TERM=NONE
BASEDIR=$(dirname "$0" | xargs -n1 realpath)
PATH_PIPELINECONFIG=${BASEDIR}/src/config/global
PATH_DATABASECONFIG=${BASEDIR}/src/config/gencodev27lift37
source $PATH_PIPELINECONFIG
source $PATH_DATABASECONFIG
BASEDIR=$(dirname "$0" | xargs -n1 realpath)
source ${BASEDIR}/input.config

#3.INPUT
##contents of input.config:
#IDT=DN20N0002ADZAA02
#IDN=DN20N0002BDZAA02
#IDR=RN20N0002RNZAA09
#OUTDIR=./
#IDPair=${IDT}-VS-${IDN}
#
#DNA_VAF=somatic/overlap/${IDPair}/${IDPair}.somatic.vaf.filter
#	#突变的频率、Ref/Alt深度信息，5列：
#	#Chr Pos VAF Ref Alt
#	#chr1    1387785 0.0221  574 13
#AAchange=somatic/overlap/${IDPair}/${IDPair}.snv.sindel.AAchange.xls
#	#突变的注释信息，尤其是氨基酸变化。至少要求12列，且为特定注释格式
#	#Chr Pos Ref Alt Impact  Gene    Effect  Transcript  Biotype cDNA    Protein Isoforms
#GERMSNP=somatic/varscan/${IDPair}/${IDPair}.snv.filt.Germline
#	#至少10列突变列表，主要关注深度，无需后面列的注释
#	#chrom   position    ref var normal_reads1   normal_reads2   normal_var_freq normal_gt   tumor_reads1    tumor_reads2
#PHASING=somatic/varscan/${IDPair}/${IDPair}.phasing.info
#	#format as following, may be empty:
#	#1       chr1    167038333,167038334     T,T
#	#2       chr10   16996392        AA
#CLONAL=somatic/clonal/${IDPair}/${IDPair}.mat.result.txt
#	#汇总了亚克隆群中各突变的频率的矩阵
#	#Chr Start   Freq    ref alt totalq  minusq  purity  p_ccf   ccf q   p_ccf_binom ccf_binom   q_binom p_prop_test q_prop_test CCF_0 ...
#HLAF=hla/overlap/${IDPair}/${IDPair}.hla1.4neo
#	#HLA分型，逗号分隔
#	#HLA-B40:01,HLA-B13:01,...
#HLALOH=somatic/hlaloh/${IDPair}/result/${IDPair}.HLALOH_summary.xls
#	#HLA_LOH杂合性缺失
#	#hla_b_40_01_02  阴性
#EXP=expression/rsem/${IDR}/${IDR}.genes.xls
#	#RNA-Seq expression output:
#	#gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM    HGNC
#	#ENSG00000000005.5_2     ENST00000373031.4_1,ENST00000485971.1_1 940.50  698.21  0.00    0.00    0.00    TNMD
#RNAVAF=mhc/expression/${IDPair}/${IDPair}.rna.vaf
#	#RNA-Seq mutation table
#	#format:chr pos AF DP.Ref DP.Alt
#	#1       90049383        0.2593  20      7
#RNAFusionBed=variant/starfusion/${IDR}/star-fusion.fusion_prediction.bedpe
#	#RNA-fusion output, from Star-fusion
#	#format:
#	#7   102665509   -1  7   -1  117407230   FBXL13>>CTTNBP2 6   -   -   0.1409

#4.Pipeline
#4.1.peptide_v2.s1.sh
OUTPUT=${OUTDIR}/mhc/peptide/${IDPair}
SAMPLE=${IDPair}
mkdir -p $OUTPUT
perl $UTIL_PEP_NMANNO $AAchange > $OUTPUT/$SAMPLE.snp.indel.txt
	#提取简化的突变信息列表，9列
perl $UTIL_PEP_MUTPEP2 $OUTPUT/$SAMPLE.snp.indel.txt $DB_REFGENE $DB_TRANS_FASTA $FAFILE $OUTPUT/peptide $GERMSNP $PHASING
	#输出peptide*系列文件。也做了一些过滤，如要求Germline AF>20%；

#4.2.peptide_v2.s2.sh: only first part of it.
INPUT=${OUTDIR}/mhc/peptide/${IDPair}
OUTPUT=${OUTDIR}/mhc/peptide/${IDPair}
SAMPLE=${IDPair}
sed 's/\*$//' $INPUT/peptide.result.txt > $OUTPUT/$SAMPLE.result.txt
$COMPARE_FILE $INPUT/peptide.result.txt $OUTPUT/$SAMPLE.result.txt
if [ -x $OUTPUT/split_MT ]; then rm -r $OUTPUT/split_MT ; fi
mkdir -p $OUTPUT/split_MT
split --numeric-suffixes=1 -l 100 -a 4 $OUTPUT/$SAMPLE.result.txt $OUTPUT/split_MT/$SAMPLE. --additional-suffix=.peptide
sed 's/\*$//' $INPUT/peptide.widetype.txt > $OUTPUT/$SAMPLE.widetype.txt
$COMPARE_FILE $INPUT/peptide.widetype.txt $OUTPUT/$SAMPLE.widetype.txt
if [ -x $OUTPUT/split_WT ]; then rm -r $OUTPUT/split_WT ; fi
mkdir -p $OUTPUT/split_WT
split --numeric-suffixes=1 -l 100 -a 4 $OUTPUT/$SAMPLE.widetype.txt $OUTPUT/split_WT/$SAMPLE. --additional-suffix=.peptide


#4.2.1. re-designed core steps, combining netMHCpan.medium.sh & TAP.medium.sh, and remove the dependence to in-house toolsets
HLAI=`cat $HLAF`
for PTYPE in WT MT ; do
	for MTPEP in `ls ${OUTDIR}/mhc/peptide/${IDPair}/split_${PTYPE}/${IDPair}*.peptide`; do
		SGE_TASK_ID=`basename $MTPEP|grep -oP '\d+.peptide'|sed 's/0//g;s/.peptide//'`
		for len in {8..11}; do
			OUTPUT=${OUTDIR}/mhc/netMHCpan/${IDPair}/$PTYPE/$len
			#netMHC
			mkdir -p $OUTPUT
			$PROG_NETMHC -a $HLAI -l $len -f $MTPEP > $OUTPUT/$SAMPLE.$SGE_TASK_ID.netMHCpan
			#TAP
			mkdir -p $OUTPUT/$SGE_TASK_ID
			perl $UTIL_TAP_FORMAT $MTPEP $OUTPUT/$SGE_TASK_ID > $OUTPUT/$SGE_TASK_ID/fa.list
			for fa in `cat $OUTPUT/$SGE_TASK_ID/fa.list`; do
				fa_name=`basename $fa | sed 's/\./_/g'`
				for hla in `sed 's/,/\n/g' $HLAF`; do
					python2 $UTIL_TAP_PREDICT --method netctlpan --allele $hla --length $len --noplot $fa > $OUTPUT/$SGE_TASK_ID/$fa_name.$hla.$len.result
				done
			done
		done
	done
done


#4.2.2.merge_netMHCpan.sh
INDIR=$OUTDIR/mhc/netMHCpan/${IDPair}
PEPID=$OUTDIR/mhc/peptide/${IDPair}/peptide.peptideID.txt
OUTPUT=${OUTDIR}/mhc/netMHCpan/${IDPair}
SAMPLE=${IDPair}
cat $INDIR/MT/*/*.netMHCpan > $OUTPUT/$SAMPLE.netMHCpan
cat $INDIR/WT/*/*.netMHCpan > $OUTPUT/$SAMPLE.widetype.netMHCpan
perl $UTIL_NETMHC_SELECT $OUTPUT/$SAMPLE.widetype.netMHCpan $PEPID $OUTPUT/$SAMPLE.netMHCpan > $OUTPUT/$SAMPLE.select.netMHCpan
perl $UTIL_NETMHC_SELECTALL $OUTPUT/$SAMPLE.widetype.netMHCpan $PEPID $OUTPUT/$SAMPLE.netMHCpan > $OUTPUT/$SAMPLE.select.netMHCpan.total
perl $UTIL_NEOANNO $PEPID $OUTPUT/$SAMPLE.select.netMHCpan > $OUTPUT/$SAMPLE.select.netMHCpan.temp
mv $OUTPUT/$SAMPLE.select.netMHCpan.temp $OUTPUT/$SAMPLE.select.netMHCpan
#rm -fr $INPUT/MT $INPUT/WT #clean

#4.3.merge_TAP.sh
INPUT=${OUTDIR}/mhc/TAP/${IDPair}
OUTPUT=${OUTDIR}/mhc/TAP/${IDPair}
SAMPLE=${IDPair}
mkdir -p $OUTPUT
find  $INPUT/WT -name '*.result' -print > $OUTPUT/wt.list
find  $INPUT/MT -name '*.result' -print > $OUTPUT/mt.list
perl $UTIL_TAP_CAT $OUTPUT/wt.list $OUTPUT/mt.list > $OUTPUT/NetCTLpan.result.txt
find  $INPUT/MT -name '*.result' -delete
find  $INPUT/WT -name '*.result' -delete
find  $INPUT/WT -name '*.fa' -delete
find  $INPUT/MT -name '*.fa' -delete


#4.4.capneo_anno.sh, neo annotation
NEOANTIGEN_RAW=${OUTDIR}/mhc/netMHCpan/${IDPair}/${IDPair}.select.netMHCpan.total
PEPTIDE_ID=${OUTDIR}/mhc/peptide/${IDPair}/peptide.peptideID.txt
DNA_VAF=${DNA_VAF}
TAP=${OUTDIR}/mhc/TAP/${IDPair}/NetCTLpan.result.txt
MUTANNO=${AAchange}
HLALOH=${HLALOH}
OUTPUT=${OUTDIR}/mhc/MHCIanno/${IDPair}
SAMPLE=${IDPair}
INFO=${OUTDIR}/mhc/peptide/${IDPair}/peptide.details
mkdir -p $OUTPUT
perl $UTIL_NEOANNO $PEPTIDE_ID $NEOANTIGEN_RAW > $OUTPUT/$SAMPLE.anno
perl $UTIL_NEO_TRANS $DB_PROTEIN $DB_TRANID_TURN $OUTPUT/$SAMPLE.anno > $OUTPUT/$SAMPLE.transcript
perl $UTIL_OVERLAP_2 <(perl $UTIL_NEOUNIQ $DB_PEPFA <(cut -f 3 $OUTPUT/$SAMPLE.transcript | sort -u) 0 | awk '$2>0 || $1=="Peptide"') $OUTPUT/$SAMPLE.transcript 2 > $OUTPUT/$SAMPLE.transcript.fp
perl $UTIL_NEO_ANCHOR $DB_ANCHORSITE $OUTPUT/$SAMPLE.transcript.fp > $OUTPUT/$SAMPLE.transcript.fp.anchor
perl $UTIL_NEO_ANNOTAP $TAP $OUTPUT/$SAMPLE.transcript.fp.anchor > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP
perl $UTIL_NEO_HOT $PATH_HOTSPOT/recurrent.position.stat.txt $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot
perl $UTIL_NEO_KNOW $DB_KNOWN_NEOEPITOPE $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope
perl $UTIL_NEO_ANNOVAF $DNA_VAF $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf
#Analysis without RNA-Seq, skipped:
#perl $UTIL_OVERLAP_2 <(cut -f 1,8 $DB_GENE_EXP) $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf 8 > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.CCLEexpStatus
#perl $UTIL_NEO_SORT $MUTANNO $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.CCLEexpStatus > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.CCLEexpStatus.sort.txt
#$UTIL_NEO_CDNA $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.CCLEexpStatus.sort.txt $INFO > $OUTPUT/$SAMPLE.HQneo.cdna.xls
#awk 'BEGIN{FS="\t";OFS="\t"}{if($39==3){$39="high";}else if($39==2){$39="medium";}else if($39==1){$39="low";}else if($39==0){$39="none";}} {if($2<7 || $1=="NeoRank"){print $0}}' $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.CCLEexpStatus.sort.txt | cut -f 1,3,6,7,9-17,20,27,38,39 |awk '$2>0' |sed 's/p\.//g' > $OUTPUT/$SAMPLE.final.report.txt
#if [ -e $HLALOH ]; then
#	$UTIL_NEO_FILTERLOH $HLALOH $OUTPUT/$SAMPLE.final.report.txt > $OUTPUT/$SAMPLE.noLOH.final.report.txt
#fi

##*.DlNeoantigen.sh: skipped
#hlatype=${HLAF}
#mhc_peptide_res=${OUTDIR}/mhc/peptide/${IDPair}/${IDPair}.result.txt
#mhc_peptide_wt=${OUTDIR}/mhc/peptide/${IDPair}/${IDPair}.widetype.txt
#mhc_peptide_id=${OUTDIR}/mhc/peptide/${IDPair}/peptide.peptideID.txt
#OUTPUT=${OUTDIR}/mhc/DlNeoantigen/${IDPair}/
#SAMPLE=${IDPair}
#mkdir -p $OUTPUT
### program path
##source $PATH_PIPELINECONFIG
#export PYTHONPATH=/mnt/nfs/software/share/conda/miniconda3/envs/python3.6/lib-oldcpu/:$PYTHONPATH
### database path
#source $PATH_DATABASECONFIG
#conda_dir=`dirname $ENV_CONDA`
#source $ENV_CONDA python3.6
#export KERAS_BACKEND=theano
#python3 $UTIL_DlNeoantigen \
#-i $hlatype \
#-m $mhc_peptide_res \
#-w $mhc_peptide_wt \
#-p $mhc_peptide_id \
#-s $SAMPLE \
#-o $OUTPUT
##source $conda_dir/deactivate
##*.DlNeoantigen_anno.sh: skipped

#4.5.FusionNeo.sh (Only analysis RNA-Seq by star-fusion)
BEDPE=${RNAFusionBed} ##input: bed file from star-fusion
HLA=${HLAF}
OUTPUT=${OUTDIR}/mhc/FusionNeo/${IDR}
mkdir -p $OUTPUT
sed 's/$/\tNA\tNA\tNA\tNA/' $HLA |sed 's/,/\tNA\tNA\tNA\tNA\n/g' > $OUTPUT/HLA.input
$PROG_FUSIONNEO/integrate-neo.py -t $OUTPUT/HLA.input -f $BEDPE -r $DB_NONCHR_FA -g $DB_GENEPRED -k -o $OUTPUT
if [ -e $OUTPUT/result.bedpe ]; then
	awk -F"\t" 'BEGIN{OFS="\t"}{if($15=="NA"){$16=$14;$15=$13;$14=$12};if($17=="NA"){$17=0};if($18=="NA"){$18=0};print}' $OUTPUT/result.bedpe>$OUTPUT/result.bedpe.tmp
	mv $OUTPUT/result.bedpe.tmp $OUTPUT/result.bedpe
fi


#4.6.RNASNV
#4.6.1.varscan_single_s1.sh #it takes RNA-Seq bam, skipped
#4.6.2.rnavaf.sh #skipped
#4.6.3.mut_rsem.sh
MUT=${AAchange}
EXP=${EXP}
VAF=${RNAVAF}
OUTPUT=${OUTDIR}/mhc/expression/${IDPair}
SAMPLE=${IDPair}
mkdir -p $OUTPUT
perl $UTIL_NEO_MUTRSEM -m $MUT -e $EXP -v $VAF -o $OUTPUT/$SAMPLE.mut.exp.xls

#4.7.Report
#4.7.1.sort_neoI_exp.sh, TruNeo ranking
NEO=${OUTDIR}/mhc/MHCIanno/${IDPair}/${IDPair}.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf
EXP=${OUTDIR}/mhc/expression/${IDPair}/${IDPair}.mut.exp.xls
CLONAL=${CLONAL}
MUTANNO=${AAchange}
HLALOH=${HLALOH}
OUTPUT=${OUTDIR}/mhc/MHCIanno/${IDPair}
SAMPLE=${IDPair}
INFO=${OUTDIR}/mhc/peptide/${IDPair}/peptide.details
mkdir -p $OUTPUT
perl $UTIL_NEO_ANNOTPM $EXP $NEO > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus
perl $UTIL_NEO_SORT $MUTANNO $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus.sort.txt
perl $UTIL_NEO_ANNOCLONAL $CLONAL $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus.sort.txt > $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus.sort.clonal.txt
$UTIL_NEO_CDNA $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus.sort.clonal.txt $INFO > $OUTPUT/$SAMPLE.HQneo.cdna.xls
awk 'BEGIN{FS="\t";OFS="\t"}{if($39==3){$39="high";}else if($39==2){$39="medium";}else if($39==1){$39="low";}else if($39==0){$39="none";}} {print $0}' $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus.sort.clonal.txt | cut -f 1,3,6,7,9-17,20,27,38- |awk '$2>0' |sed 's/p\.//g' > $OUTPUT/$SAMPLE.final.report.txt
if [ -e $HLALOH ]; then
	$UTIL_NEO_FILTERLOH $HLALOH $OUTPUT/$SAMPLE.final.report.txt > $OUTPUT/$SAMPLE.noLOH.final.report.txt
fi
perl $UTIL_NEO_FORMAT1 -a $OUTPUT/$SAMPLE.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf.RSEMexpStatus.sort.clonal.txt -t RNA -o1 $OUTPUT/$SAMPLE.MHCI.summary.neoantigen.xls -o2 $OUTPUT/$SAMPLE.MHCI.summary.neoantigen.WBSB.xls
perl $UTIL_NEO_FORMAT2 -i1 $OUTPUT/$SAMPLE.MHCI.summary.neoantigen.xls -i2 $OUTPUT/$SAMPLE.MHCI.summary.neoantigen.WBSB.xls -t RNA -o1 $OUTPUT/$SAMPLE.MHCI.summary.xls -o2 $OUTPUT/$SAMPLE.MHCI.summary.WBSB.xls


#4.7.2./mnt/nfs/pipeline/yucecloud/modular/template/report/neoantigen.v2.sh
#Result filter
MUTNEO=${OUTDIR}/mhc/MHCIanno/${IDPair}/${IDPair}.noLOH.final.report.txt
##${INPUT[1]}=mhc/FusionNeo/${IDPair}/result.bedpe #skipped, form DNA
FFusion=${OUTDIR}/mhc/FusionNeo/${IDR}/result.bedpe
OUTPUT=${OUTDIR}/report/neoantigen/${IDPair}
SAMPLE=${IDPair}
#PARAM=${OUTDIR}/mhc/DlNeoantigenAnno/${IDPair}/${IDPair}.noLOH.final.report.txt #skipped
INFO=MHCI.final.report.txt
mkdir -p $OUTPUT
if [ -e ${FFusion} ]; then
  $UTIL_NEO_MERGE ${MUTNEO} ${FFusion} > $OUTPUT/$SAMPLE.$INFO
  awk -F"\t" 'NR == 1; NR > 1 {print $0 | "sort -rnk2"}' $OUTPUT/$SAMPLE.$INFO|awk -F"\t" 'BEGIN{a=0;OFS="\t"}{if(/^NeoRank/){print;next}$1="";print(++a$0)}' - >$OUTPUT/$SAMPLE.$INFO.tmp
  mv $OUTPUT/$SAMPLE.$INFO.tmp $OUTPUT/$SAMPLE.$INFO
else
  cp $MUTNEO $OUTPUT/$SAMPLE.$INFO
fi
if [ -f $OUTPUT/MUTNEO.tmp ]; then
	rm $OUTPUT/MUTNEO.tmp
fi
#限定config文件

#4.7.3./mnt/nfs/pipeline/yucecloud/modular/template/mhc/merge2xlsx_v2.sh #skipped
#generate spreadsheet report

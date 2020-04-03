#!/bin/bash
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
# Description: This is the main scripts of TruNeo, which 
# aimed at predicting and prioritizing neoantigens from 
# DNA Variation & RNA-Seq.
#
#=============================================

#1.Prerequisite
# see details in README.txt

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
##Please configure the content of input.config:

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


#4.7.2.neoantigen.v2.sh
#Result filter
MUTNEO=${OUTDIR}/mhc/MHCIanno/${IDPair}/${IDPair}.noLOH.final.report.txt
##${INPUT[1]}=mhc/FusionNeo/${IDPair}/result.bedpe #skipped, form DNA
FFusion=${OUTDIR}/mhc/FusionNeo/${IDR}/result.bedpe
OUTPUT=${OUTDIR}/report/neoantigen/${IDPair}
SAMPLE=${IDPair}
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

#The final result is at $OUTDIR/report/neoantigen/

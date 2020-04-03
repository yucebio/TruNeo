#!/bin/bash
export R_LIBS=/mnt/nfs/database/global/yucemed/lohhla/env/lib/R:$R_LIBS
export PATH="/mnt/nfs/database/global/yucemed/lohhla/env/anaconda3/bin:/mnt/nfs/software/share/jre1.8.0_91/bin:/mnt/nfs/software/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/mnt/nfs/database/global/yucemed/lohhla/env/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/opt/dell/srvadmin/bin:$PATH"
export LD_LIBRARY_PATH="/mnt/nfs/database/global/yucemed/lohhla/env/lib:/mnt/nfs/software/lib:/opt/MonitorSoftware/lib:/usr/sfw/lib:/usr/local/lib:$LD_LIBRARY_PATH"
#module load BEDTools/2.26.0-foss-2016b
#module load SAMtools/1.3.1-foss-2016b
#module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3 # for R package requirements
#module load novoalign/3.07.00
#module load TracerX-Picard-GATK/0.1-Java-1.7.0_80
#module load Jellyfish/2.2.6-foss-2016b
# If these commands are not already in your path, they must be added or pointed to
#alias samtools=/path/to/samtools
#alias jellyfish=/mnt/nfs/user/renqian/bin/jellyfish/bin/jellyfish
#alias novoindex=/mnt/nfs/user/renqian/bin/novocraft/novoindex
#alias bedtools=/path/to/bedtools
sample=$1
outtput=$2
normal=$3
bam=$4
hlas=$5
solution=$6
Rscript /mnt/nfs/database/global/yucemed/lohhla/LOHHLAscript.R \
        --patientId $sample \
        --outputDir $outtput \
        --normalBAMfile $normal \
        --BAMDir $bam \
        --hlaPath $hlas \
        --HLAfastaLoc /mnt/nfs/database/global/yucemed/lohhla/data/hla_all.fasta \
        --CopyNumLoc $solution \
        --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE \
        --HLAexonLoc /mnt/nfs/database/global/yucemed/lohhla/data/hla.dat \
        --gatkDir /mnt/nfs/database/global/yucemed/lohhla/picard-tools-1.119/ \
        --novoDir /mnt/nfs/software/share/novocraft_v3.03/

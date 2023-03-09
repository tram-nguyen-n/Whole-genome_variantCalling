#!/bin/bash -l
#SBATCH --job-name=fsjwgsP3
#SBATCH --output=fsjwgsP3-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32000
#SBATCH --time=16:00:00
#SBATCH --partition=regular,long7,long30
#SBATCH --exclude=cbsubscb02
#SBATCH --account=bscb02
#SBATCH --array=1-232

# phase3: GATK HaplotypeCaller
# # current unique sample count: 232
# exclude the longer samples: 82,86,208,228

# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

allSamples=( `cat full_unique_sample_list.txt` )
mySample=${allSamples[${SLURM_ARRAY_TASK_ID}-1]}

allSex=( `cat full_unique_sample_sex_list.txt` )
mySex=${allSex[${SLURM_ARRAY_TASK_ID}-1]}

myOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-10-01/per_sample_output/${mySample}

echo $HOSTNAME
echo sample name
echo $mySample
echo output directory
echo $myOutDir
echo task index
echo $SLURM_ARRAY_TASK_ID
echo $newdir

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir

echo current directory /workdir/$USER/$newdir

## create tmp directory for java
mkdir tmp

## copy required files
/programs/bin/labutils/mount_server cbsuclarkfs1 /storage
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/per_sample_output/${mySample}/${mySample}_phase2.bam .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/FSJv2_reference/*.* .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/contigs.txt .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/contigs_auto.txt .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/contigs_Z.txt .

pwd


## edit bam file scaffold names and index edited bam
#samtools view -H ${mySample}_phase2.bam | sed -e 's/;HRSCAF=[0-9]\+//' | samtools reheader - ${mySample}_phase2.bam > ${mySample}.bam
# SEPT 15, 2020 -- no longer need to edit header line for bam since now mapped to edited reference from the beginning
mv ${mySample}_phase2.bam ${mySample}.bam
samtools index ${mySample}.bam ${mySample}.bai

## HaplotypeCaller
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH


echo
echo start HaplotypeCaller
date
if [ "$mySex" == "M" ]; then
  echo skipping male sample
  echo running male haplotypecaller
  for myContig in `cat contigs.txt`; do
    java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -Xmx32g -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R FSJ.v2_1.fasta -I ${mySample}.bam -L ${myContig} -nct 8 -secondsBetweenProgressUpdates 600 --emitRefConfidence GVCF -o ${mySample}.${myContig}.g.vcf
  done
else
  echo running female autosome haplotypecaller
  for myContig in `cat contigs_auto.txt`; do
    java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -Xmx32g -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R FSJ.v2_1.fasta -I ${mySample}.bam -L ${myContig} -nct 8 -secondsBetweenProgressUpdates 600 --emitRefConfidence GVCF -o ${mySample}.${myContig}.g.vcf
  done
  echo running female Z haplotypecaller
  myContig=`cat contigs_Z.txt`
  java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -Xmx32g -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R FSJ.v2_1.fasta -I ${mySample}.bam -L ${myContig} -ploidy 1 -nct 8 -secondsBetweenProgressUpdates 600 -ERC GVCF -o ${mySample}.${myContig}.g.vcf
fi
echo completed HaplotypeCaller
date



## copy output to storage
mkdir -p $myOutDir/per_contig_gvcf
cp ${mySample}.*.g.vcf* $myOutDir/per_contig_gvcf


## clean up
cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

#!/bin/bash -l
#SBATCH --job-name=fsjwgsP2
#SBATCH --output=fsjwgsP2-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48000
#SBATCH --time=24:00:00
#SBATCH --partition=regular,long7,long30
#SBATCH --exclude=cbsubscb02,cbsubscb17,cbsubscb16
#SBATCH --account=bscb02
#SBATCH --array=1-291

# phase2: post-alignment clean-up
# run various Picard and GATK tools
# run per sample, combining files from multiple runs in MarkDuplicates step
# full_unique_sample_list.txt is a single column with sample names
# current unique sample count: 233-291 (final batch of samples)
# omitting base recalibration step

# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

allSamples=( `cat full_unique_sample_list_2022.txt` )
mySample=${allSamples[${SLURM_ARRAY_TASK_ID}-1]}

myOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2022-Jan/per_sample_output/${mySample}
myOldOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/per_sample_output/${mySample}

echo $HOSTNAME
echo sample name
echo $mySample
echo output directory
echo $myOutDir
echo task index
echo $SLURM_ARRAY_TASK_ID
echo $newdir

mkdir -p /workdir/$USER/$newdir

# copy run info
cp run_info_2022.txt /workdir/$USER/$newdir

cd /workdir/$USER/$newdir

## create tmp directory for java
mkdir tmp

## copy required files
/programs/bin/labutils/mount_server cbsuclarkfs1 /storage
cp $myOutDir/*/${mySample}_*.sorted.bam .
cp $myOldOutDir/*/${mySample}_*.sorted.bam . ##grab older samples if they exist for reruns

pwd

echo
ls *.sorted.bam
echo

## identify runs for current sample
grep "${mySample}[[:space:]]" run_info_2022.txt > c_info.txt
runCount=( `grep -c "${mySample}[[:space:]]" run_info_2022.txt` )


## AddOrReplaceReadGroups
# NOTE: if RGID set properly in phase1, this step can be skipped!
for runIdx in $(seq 1 ${runCount}); do
  myRun=`sed -n "${runIdx}p" c_info.txt | cut -f 2`
  myLane=`sed -n "${runIdx}p" c_info.txt | cut -f 3`
  runOnly=`echo $myRun | sed 's/_.*//'`
  echo
  echo Output file name is=${mySample}_${myRun}_sorted.bam
  echo RGID=${runOnly}.${myLane}
  echo RGLB=${mySample}
  echo RGPL=illumina RGPU=${runOnly}.${myLane}.${mySample}
  echo RGSM=${mySample}
  echo
  echo start AddOrReplaceReadGroups for ${myRun}
  date
  java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/picard-tools-2.8.2/picard.jar AddOrReplaceReadGroups I=${mySample}_${myRun}.sorted.bam O=${mySample}_${myRun}_sorted.bam RGID=${runOnly}.${myLane} RGLB=${mySample} RGPL=illumina RGPU=${runOnly}.${myLane}.${mySample} RGSM=${mySample} SORT_ORDER=coordinate CREATE_INDEX=true QUIET=true VERBOSITY=WARNING
  echo completed AddOrReplaceReadGroups for ${myRun}
  date
done


## MarkDuplicates
# During this step -- we will also combine all .bam files produced for each sample (because some samples were run over multiple lanes and have different RGID)
# GATK: run the initial steps of the pre-processing workflow (mapping and sorting) separately on each individual read group. Then we merge the data to produce a single BAM file for each sample (aggregation); this is done conveniently at the same time that we do the duplicate marking, by running Mark Duplicates on all read group BAM files for a sample at the same time.
# create array variable with all input files
ls *_sorted.bam > bam_list.txt
awk '{print "INPUT="$1}' bam_list.txt > input_bam_list.txt
readarray -t bam_array < input_bam_list.txt
echo ${bam_array[@]}
# run
echo
echo start MarkDuplicates and consolidating sample bam files
date
java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates $(echo ${bam_array[@]}) OUTPUT=${mySample}_markDup.bam METRICS_FILE=${mySample}_markDup_metrics.txt CREATE_INDEX=true ASSUME_SORTED=true QUIET=true VERBOSITY=WARNING
echo completed MarkDuplicates and consolidating sample bam files
date


## intermediate clean up
rm *sorted.bam


## FixMate
echo
echo start FixMate
date
java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -Xmx48g -jar /programs/picard-tools-1.109/FixMateInformation.jar INPUT=${mySample}_markDup.bam OUTPUT=${mySample}_phase2.bam SO=coordinate CREATE_INDEX=true QUIET=true VERBOSITY=WARNING
echo completed FixMate
date


## run qualimap
echo
echo start qualimap
date
/programs/qualimap_v2.2.1/qualimap bamqc -bam ${mySample}_phase2.bam -c -nw 400 -hm 3 --java-mem-size=4G
echo completed qualimap
date


## list all files
ls -l


## copy output to storage
cp ${mySample}_phase2.bam $myOutDir
cp -r ${mySample}_phase2_stats $myOutDir
cp ${mySample}_markDup_metrics.txt $myOutDir


## clean up
cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

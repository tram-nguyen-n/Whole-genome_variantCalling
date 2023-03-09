#!/bin/bash -l
#SBATCH --job-name=fsjwgsP1
#SBATCH --output=fsjwgsP1-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --partition=regular,long7,long30
#SBATCH --exclude=cbsubscb02
#SBATCH --account=bscb02
#SBATCH --array=1-357

# phase1: pre-processing FSJ WGS
# For us, 357 represents the number of all unique sequence runs that we have
# run fastqc, trimmomatic, BWA, and sort output bam

# get your array number from the fsj_wgs_sample_list.txt
# may need to exclude=cbsubscb02
# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

mySample=`sed -n "${SLURM_ARRAY_TASK_ID}p" fsj_wgs_sample_list.txt | cut -f 1`
myRun=`sed -n "${SLURM_ARRAY_TASK_ID}p" fsj_wgs_sample_list.txt | cut -f 2`
myOutStr=${mySample}_${myRun}
myOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/per_sample_output/${mySample}/${myRun}

echo $HOSTNAME
echo sample name
echo $mySample
echo run ID
echo $myRun
echo output directory
echo $myOutDir
echo output file prefix
echo $myOutStr
echo task index
echo $SLURM_ARRAY_TASK_ID
echo $newdir

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir

/programs/bin/labutils/mount_server cbsuclarkfs1 /storage
cp /fs/cbsuclarkfs1/storage/ejc87/fsj/tram_wgs/rawdata/${mySample}_C*${myRun}*.fq.gz .
cp /fs/cbsuclarkfs1/storage/ejc87/fsj/tram_wgs/pre-processing/2020-01-28/NEBNext_adapters.fa . #file of adapter sequences to remove
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/FSJv2_reference/*.* .

mkdir -p $myOutDir

pwd

## fastqc
echo running fastqc
export PATH=/programs/FastQC-0.11.8:$PATH
fastqc --quiet *_1.fq.gz
fastqc --quiet *_2.fq.gz
mv *_fastqc.* $myOutDir

## trimmomatic
echo running trimmomatic
java -jar /programs/trimmomatic/trimmomatic-0.36.jar PE -threads ${SLURM_NTASKS} ./*_1.fq.gz ./*_2.fq.gz -baseout ${myOutStr}.fq.gz ILLUMINACLIP:NEBNext_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10
echo finish trimmomatic


## alignment
echo running alignment BWA MEM
/programs/bwa-0.7.13/bwa mem -M -t ${SLURM_NTASKS} FSJ.v2_1.fasta ./${myOutStr}_1P.fq.gz ./${myOutStr}_2P.fq.gz > ./${myOutStr}.sam


## get alignment stats from BWA
samtools flagstat ./${myOutStr}.sam > ./${myOutStr}.log # currently on the local directory on node your task sent to (deleted after your job!)
mv *.log $myOutDir


## sam to bam
samtools view -bS ${myOutStr}.sam > ${myOutStr}.bam


## sort bam
samtools sort ${myOutStr}.bam -o ${myOutStr}.sorted.bam


## copy output to storage
mv *.sorted.bam $myOutDir

cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

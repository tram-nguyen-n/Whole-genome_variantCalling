#!/bin/bash -l
#SBATCH --job-name=fsjwgsP4
#SBATCH --output=fsjwgsP4-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16000
#SBATCH --time=24:00:00
#SBATCH --partition=regular,long7,long30
#SBATCH --exclude=cbsubscb02
#SBATCH --account=bscb02
#SBATCH --array=1-866

# phase4: GATK GenotypeGVCFs
# DO NOT include Z: contig 30
# contig count: all shorter contigs NOT 1-5,11-13,629,866
# here our array number represents the number of contigs we want to run 

# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

allContigs=( `cat contigs.txt` )
myContig=${allContigs[${SLURM_ARRAY_TASK_ID}-1]}

allSex=( `cat full_unique_sample_sex_list.txt` )
mySex=${allSex[${SLURM_ARRAY_TASK_ID}-1]}

#MAKE SURE TO CREATE THIS OUTPUT DIRECTORY ALREADY
myOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-10-01/per_contig_vcfs

echo $HOSTNAME
echo contig name
echo $myContig
echo task index
echo $SLURM_ARRAY_TASK_ID
echo $newdir

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir

## create tmp directory for java
mkdir tmp

## copy required files
/programs/bin/labutils/mount_server cbsuclarkfs1 /storage
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/per_sample_output/*/per_contig_gvcf/*.${myContig}.g.vcf* .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/FSJv2_reference/*.* .
cp ~/fsj_wgs20x/2021-04-14/full_unique_sample_list.txt .


## get input file string
awk -v myVar="$myContig" '{print "--variant "$1"."myVar".g.vcf"}' full_unique_sample_list.txt > input_gvcf_list.txt
readarray -t gvcf_array < input_gvcf_list.txt
echo ${gvcf_array[@]}


## GenotypeGVCFs
echo
echo start GenotypeGVCFs
date
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH
java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -Xmx16g -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -L ${myContig} -secondsBetweenProgressUpdates 600 -R FSJ.v2_1.fasta $(echo ${gvcf_array[@]}) -o unfiltered_variants_fsj_wgs.${myContig}.vcf
echo completed GenotypeGVCFs
date


## copy output to storage
cp unfiltered_variants_fsj_wgs.${myContig}.vcf* $myOutDir


## clean up
cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

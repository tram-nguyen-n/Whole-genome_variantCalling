#!/bin/bash -l
#SBATCH --job-name=fsjwgsP7-allsites
#SBATCH --output=fsjwgsP7-allsites-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16000
#SBATCH --time=4:00:00
#SBATCH --partition=short,regular
#SBATCH --account=bscb02
#SBATCH --array=1-878

# phase6: Vcftools filter invariants -- quick contigs
# contig count: 1-878

# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

allContigs=( `cat contigs.txt` )
myContig=${allContigs[${SLURM_ARRAY_TASK_ID}-1]}

myOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/per_contig_filtered_invariant_vcfs2

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
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/per_contig_invariant_vcfs2/invariant_${myContig}.vcf.gz* .
#cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/FSJv2_reference/*.* .

pwd

## Filter with vcftools
echo
echo start VCFtools filtering
date

## additional filtering checks -- based on variants depth
echo filtering 1 based on variant filters
vcftools --gzvcf invariant_${myContig}.vcf.gz --minGQ 20 --minQ 30 --recode --remove-filtered-all --out temp1.${myContig}
vcftools --vcf temp1.${myContig}.recode.vcf --min-meanDP 10 --max-meanDP 121 --recode --remove-filtered-all --out temp2.${myContig}
vcftools --vcf temp2.${myContig}.recode.vcf --max-missing 0.85 --recode --remove-filtered-all --out filter1.allsites.${myContig}
echo additional filter 1 count
#echo `grep -v '#' filter1.allsites.${myContig}.recode.vcf | wc -l`
echo


## based on invariant stats
echo filtering 2 based on chrom 1 stats
vcftools --gzvcf invariant_${myContig}.vcf.gz --minGQ 20 --minQ 30 --recode --remove-filtered-all --out temp1.${myContig}
vcftools --vcf temp1.${myContig}.recode.vcf --min-meanDP 10 --max-meanDP 55 --recode --remove-filtered-all --out temp2.${myContig}
vcftools --vcf temp2.${myContig}.recode.vcf --max-missing 0.85 --recode --remove-filtered-all --out filter2.allsites.${myContig}
echo additional filter 2 count
#echo `grep -v '#' filter2.allsites.${myContig}.recode.vcf | wc -l`
echo

## low depth
echo filtering 3 low depth
vcftools --gzvcf invariant_${myContig}.vcf.gz --minGQ 20 --minQ 30 --recode --remove-filtered-all --out temp1.${myContig}
vcftools --vcf temp1.${myContig}.recode.vcf --min-meanDP 10 --max-meanDP 200 --recode --remove-filtered-all --out temp2.${myContig}
vcftools --vcf temp2.${myContig}.recode.vcf --max-missing 0.85 --recode --remove-filtered-all --out filter3.allsites.${myContig}
echo additional filter 3 count
#echo `grep -v '#' filter3.allsites.${myContig}.recode.vcf | wc -l`
echo


## low depth + missingness
echo filtering 4 low missing and depth
vcftools --gzvcf invariant_${myContig}.vcf.gz --minGQ 20 --minQ 30 --recode --remove-filtered-all --out temp1.${myContig}
vcftools --vcf temp1.${myContig}.recode.vcf --min-meanDP 10 --max-meanDP 200 --recode --remove-filtered-all --out temp2.${myContig}
vcftools --vcf temp2.${myContig}.recode.vcf --max-missing 0.80 --recode --remove-filtered-all --out filter4.allsites.${myContig}
echo additional filter 4 count
#echo `grep -v '#' filter4.allsites.${myContig}.recode.vcf | wc -l`
echo


echo bgzipping all vcfs
date

bgzip -c filter1.allsites.${myContig}.recode.vcf > filter1.invariant.${myContig}.vcf.gz
bgzip -c filter2.allsites.${myContig}.recode.vcf > filter2.invariant.${myContig}.vcf.gz
bgzip -c filter3.allsites.${myContig}.recode.vcf > filter3.invariant.${myContig}.vcf.gz
bgzip -c filter4.allsites.${myContig}.recode.vcf > filter4.invariant.${myContig}.vcf.gz

echo finished zipping
date

# index them all
tabix filter1.invariant.${myContig}.vcf.gz
tabix filter2.invariant.${myContig}.vcf.gz
tabix filter3.invariant.${myContig}.vcf.gz
tabix filter4.invariant.${myContig}.vcf.gz

echo
echo copying files
echo

## copy output to storage
cp filter1.invariant.${myContig}.vcf* $myOutDir
cp filter1.invariant.${myContig}.vcf* $myOutDir
cp filter1.invariant.${myContig}.vcf* $myOutDir
cp filter1.invariant.${myContig}.vcf* $myOutDir


## clean up
cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

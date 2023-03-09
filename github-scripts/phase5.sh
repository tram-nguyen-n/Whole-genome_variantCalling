#!/bin/bash -l
#SBATCH --job-name=fsjwgsP5
#SBATCH --output=fsjwgsP5-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16000
#SBATCH --time=4:00:00
#SBATCH --partition=short,regular
#SBATCH --exclude=cbsubscb02,cbsubscb17
#SBATCH --account=bscb02
#SBATCH --array=1-878


# phase5: filter variants
# contig count: 878

# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

allContigs=( `cat contigs.txt` )
myContig=${allContigs[${SLURM_ARRAY_TASK_ID}-1]}

myOutDir=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-10-01/per_contig_GATKfiltered_vcfs

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
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-10-01/per_contig_vcfs/unfiltered_variants_fsj_wgs.${myContig}.vcf* .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/FSJv2_reference/*.* .

pwd


## SelectVariants
echo
echo start SelectVariants
date
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH
#java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T SelectVariants -L ${myContig} -secondsBetweenProgressUpdates 600 -R FSJ.v2_1.fasta -V unfiltered_variants_fsj_wgs.${myContig}.vcf -selectType SNP -o unfiltered_snps.${myContig}.vcf

java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T SelectVariants -L ${myContig} -secondsBetweenProgressUpdates 600 -R FSJ.v2_1.fasta -V unfiltered_variants_fsj_wgs.${myContig}.vcf -selectType INDEL -o unfiltered_indels.${myContig}.vcf

echo completed indels

echo start selecting All Variants

java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T SelectVariants -L ${myContig} -secondsBetweenProgressUpdates 600 -R FSJ.v2_1.fasta -V unfiltered_variants_fsj_wgs.${myContig}.vcf -o unfiltered_AllVariants.${myContig}.vcf
echo completed All Variants


echo completed SelectVariants
date


## VariantFiltration
#echo
#echo start VariantFiltration
#date
#export JAVA_HOME=/usr/local/jdk1.8.0_121
#export PATH=$JAVA_HOME/bin:$PATH

#java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T VariantFiltration \
#-L ${myContig} \
#-secondsBetweenProgressUpdates 600 \
#-R FSJ.v2_1.fasta \
#-V unfiltered_snps.${myContig}.vcf \
#--filterExpression "QD < 2.0" --filterName "filtQD" \
#--filterExpression "FS > 60.0" --filterName "filtFS" \
#--filterExpression "MQ < 40.0" --filterName "filtMQ" \
#--filterExpression "ReadPosRankSum < -8.0" --filterName "filtRP" \
#-o filtered_snps_sites.${myContig}.vcf

#echo completed VariantFiltration
#date
#echo filtered_snps_sites.${myContig}.vcf variant count
#echo `grep -v '#' filtered_snps_sites.${myContig}.vcf | wc -l`


## VariantFiltration
echo
echo start IndelFilter
date
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH

java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T VariantFiltration \
-L ${myContig} \
-secondsBetweenProgressUpdates 600 \
-R FSJ.v2_1.fasta \
-V unfiltered_indels.${myContig}.vcf \
--filterExpression "QD < 2.0" --filterName "filtQD" \
--filterExpression "FS > 60.0" --filterName "filtFS" \
--filterExpression "MQ < 40.0" --filterName "filtMQ" \
--filterExpression "ReadPosRankSum < -8.0" --filterName "filtRP" \
-o filtered_indels.${myContig}.vcf

echo completed IndelFilter
date
echo filtered_indels.${myContig}.vcf variant count
echo `grep -v '#' filtered_indels.${myContig}.vcf | wc -l`



## VariantFiltration
echo
echo start All VariantFiltration
date
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH

java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T VariantFiltration \
-L ${myContig} \
-secondsBetweenProgressUpdates 600 \
-R FSJ.v2_1.fasta \
-V unfiltered_AllVariants.${myContig}.vcf \
--filterExpression "QD < 2.0" --filterName "filtQD" \
--filterExpression "FS > 60.0" --filterName "filtFS" \
--filterExpression "MQ < 40.0" --filterName "filtMQ" \
--filterExpression "ReadPosRankSum < -8.0" --filterName "filtRP" \
-o filtered_AllVariants.${myContig}.vcf

echo completed All VariantFiltration
date
echo filtered_AllVariants.${myContig}.vcf variant count
echo `grep -v '#' filtered_AllVariants.${myContig}.vcf | wc -l`



## get FILTER field counts
#echo
#echo FILTER counts
#echo filtQD:`grep -c 'filtQD' filtered_snps.${myContig}.vcf`
#echo filtFS:`grep -c 'filtFS' filtered_snps.${myContig}.vcf`
#echo filtMQ:`grep -c 'filtMQ' filtered_snps.${myContig}.vcf`
#echo filtRP:`grep -c 'filtRP' filtered_snps.${myContig}.vcf`
#echo PASS:`grep -c 'PASS' filtered_snps.${myContig}.vcf`
#echo


echo
echo FILTER INDEL counts
echo filtQD:`grep -c 'filtQD' filtered_indels.${myContig}.vcf`
echo filtFS:`grep -c 'filtFS' filtered_indels.${myContig}.vcf`
echo filtMQ:`grep -c 'filtMQ' filtered_indels.${myContig}.vcf`
echo filtRP:`grep -c 'filtRP' filtered_indels.${myContig}.vcf`
echo PASS:`grep -c 'PASS' filtered_indels.${myContig}.vcf`
echo


echo FILTER ALL counts
echo filtQD:`grep -c 'filtQD' filtered_AllVariants.${myContig}.vcf`
echo filtFS:`grep -c 'filtFS' filtered_AllVariants.${myContig}.vcf`
echo filtMQ:`grep -c 'filtMQ' filtered_AllVariants.${myContig}.vcf`
echo filtRP:`grep -c 'filtRP' filtered_AllVariants.${myContig}.vcf`
echo PASS:`grep -c 'PASS' filtered_AllVariants.${myContig}.vcf`
echo


## additional filtering checks
#vcftools --vcf filtered_snps.${myContig}.vcf --max-missing 0.8 --min-meanDP 2 --mac 3 --max-meanDP 200 --recode --remove-filtered-all --out filter1.${myContig}.vcf
#echo additional filter 1 count
#echo `grep -v '#' filter1.${myContig}.vcf.recode.vcf | wc -l`
#echo
#vcftools --vcf filtered_snps.${myContig}.vcf --max-missing 0.90 --min-meanDP 2 --mac 3 --max-meanDP 200 --recode --remove-filtered-all --out filter2.${myContig}.vcf
#echo additional filter 2 count
#echo `grep -v '#' filter2.${myContig}.vcf.recode.vcf | wc -l`
#echo
#vcftools --vcf filtered_snps.${myContig}.vcf --max-missing 0.95 --min-meanDP 2 --mac 3 --max-meanDP 1000 --recode --remove-filtered-all --out filter3.${myContig}.vcf
#echo additional filter 3 count
#echo `grep -v '#' filter3.${myContig}.vcf.recode.vcf | wc -l`
#echo
#vcftools --vcf filtered_snps.${myContig}.vcf --max-missing 0.85 --min-meanDP 2 --max-meanDP 200 --recode --remove-filtered-all --out filter4.${myContig}.vcf
#echo additional filter 4 count
#echo `grep -v '#' filter4.${myContig}.vcf.recode.vcf | wc -l`
#echo
#max-meanDP Loci that have high mean depth are indicative of either paralogs or multicopy loci.


## REMOVE the variants that aren't passing
echo
echo start Remove variants failing GATK filtering
date
#export JAVA_HOME=/usr/local/jdk1.8.0_121
#export PATH=$JAVA_HOME/bin:$PATH

#java -Djava.io.tmpdir=/workdir/$USER/$newdir/tmp -jar /programs/bin/GATK/GenomeAnalysisTK.jar -T VariantFiltration \
#-L ${myContig} \
#-secondsBetweenProgressUpdates 600 \
#-R FSJ.v2_1.fasta \
#-V filtered_snps.${myContig}.vcf \
#--exclude-filtered \
#-o filtered_snps.${myContig}.vcf

#echo filtered_snps.${myContig}.vcf variant count
#echo `grep -v '#' filtered_snps.${myContig}.vcf | wc -l`

#awk '/^#/||$7=="PASS"' filtered_snps.${myContig}.vcf > passed_snps_only.${myContig}.vcf

awk '/^#/||$7=="PASS"' filtered_indels.${myContig}.vcf > passed_indels.${myContig}.vcf

awk '/^#/||$7=="PASS"' filtered_AllVariants.${myContig}.vcf > passed_AllVariants.${myContig}.vcf

echo completed Remove variants failing GATK filtering
date

#echo passed_snps_only.${myContig}.vcf variant count
#echo `grep -v '#' passed_snps_only.${myContig}.vcf | wc -l`

echo passed_indels.${myContig}.vcf variant count
echo `grep -v '#' passed_indels.${myContig}.vcf | wc -l`
echo
echo passed_AllVariants.${myContig}.vcf variant count
echo `grep -v '#' passed_AllVariants.${myContig}.vcf | wc -l`


## copy output to storage
#cp passed_snps_only.${myContig}.vcf* $myOutDir
cp passed_indels.${myContig}.vcf* $myOutDir/indels
cp passed_AllVariants.${myContig}.vcf* $myOutDir/all_variants


## clean up
cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

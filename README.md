# FSJ Variant Calling and Filtering
Scalable pipeline for variant discovery and quality filtering for paired-end, short-read sequences. This pipeline was developed in collaboration with Dr. Elissa Cosgrove at Cornell University. For questions regarding this pipeline, feel free to contact Tram Nguyen at tn337@cornell.edu

Here, we present the pipeline using Florida Scrub-Jay whole-genome sequences. However, this workflow can be adapted to other genomic short-read datasets. This workflow follows the GATK Best Practices protocol and is executed on Cornell's HPC clusters with a SLURM scheduler. Broadly, the pipeline consists of six phases: from alignment to a reference genome, to variant calling using GATK HaplotypeCaller, to additional quality filtering (site depth, missingness, genotyping quality, etc). Note that additional filtering thresholds are specific to our dataset and should be modified accordingly. 

# Dependencies
### Software
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.8
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.36
- [BWA](https://bio-bwa.sourceforge.net/) v0.7.13
- [Samtools](http://www.htslib.org/doc/samtools) v1.17
   - [flagstat](http://www.htslib.org/doc/samtools-flagstat.html)
   - [view](http://www.htslib.org/doc/samtools-view.html)
   - [sort](https://www.htslib.org/doc/samtools-sort.html)
- [Picard](https://broadinstitute.github.io/picard/) v2.8.2
   - [AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-)
   - [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
   - [FixMate](https://gatk.broadinstitute.org/hc/en-us/articles/360036713471-FixMateInformation-Picard-) v1.109
   - [GatherVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360037422071-GatherVcfs-Picard-)
- [QualiMap](http://qualimap.conesalab.org/) v2.2.1 
- [GATK](https://gatk.broadinstitute.org/hc/en-us) v3.8
   - [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants)
   - [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)
- [VCFtools](https://vcftools.sourceforge.net/) v0.1.16
- [BCFtools](https://samtools.github.io/bcftools/bcftools.html) v1.15.1.

# Phase 1
Phase 1: Run fastqc, trimmomatic, BWA, and sort output bam. 

Input required for this phase:
  1. Raw short read sequences with the follow file suffix, for read 1 and read 2. ```./*_1.fq.gz ./*_2.fq.gz```
  2. A tab-delimited text file with Sample Name and Run Flow Cell ID information. Below is an example:
  ``` 
SAMPLE1 HNLMNDSXX
SAMPLE2 HNK7FDSXX
SAMPLE3 HNLJLDSXX
SAMPLE4 HNLMNDSXX
  ```
  Below are helpful guides to understanding the .fq.gz file and how to extract run, flowcell, sequencer information.
  
  (https://uofabioinformaticshub.github.io/Intro-NGS-fib/notes/3-raw_data.html)
  
  (https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211)
  
  4. A file containing adapter sequences to remove in Trimmomatic.
  5. And of course, your reference genome to align to.
  
# Phase 2
Phase 2: This phase runs various Picard and GATK tools to assign read groups information to samples, mark potential PCR duplicates and if a sample was sequenced several times, this step also simultaneously combines those runs into one .bam file. Because we are working with paired-end sequences, this step also verifies mate-pair information between mates and fixes mismatches. Finally, we run Qualimap as another intermediate QC step to visualize coverage, GC content, percent of reads map, etc.

Input required for this phase:
   1. Sorted .bams from Phase 1
   2. Tab-delimited text file containing Sample ID, FlowCell ID, Sequencing Lane, and Run number. An example file is given below:
   
  ```
SAMPLE1  HNLMNDSXX       L1      1
SAMPLE2  HNK7FDSXX       L2      1
SAMPLE3  HNLJLDSXX       L1      1
SAMPLE4  HNLMNDSXX       L1      1
SAMPLE5  HNLMNDSXX       L1      1
SAMPLE6  HNK7FDSXX       L1      1
SAMPLE7  HNK7FDSXX       L1      1
SAMPLE8  HNLJLDSXX       L1      1
```

# Phase 3
Phase 3: Variant calling occurs during this step using GATK's HaplotypeCaller. To optimize efficiency, we run this phase PER CONTIG (or chromosomes in cases where you have chromosome level assemblies) and parallelize per contig runs as individual SLURM jobs. As a note, HaplotypeCaller assumes default diploid ploidy. Therefore, we separate our samples by male and females and run HaplotypeCaller with option ```-ploidy 1``` for females with one Z chromosomes and one W chromosome. If running a human dataset, this same option should be applied on the Y chromosome for males. 

Input required for this phase:
   1. Per-sample .bam files from Phase 2
   2. Your reference genome
   3. A file listing all contigs (separate out the autosomal and sex chromosomes to apply different --ploidy HaplotypeCaller options)
   4. A list of the sample names

The outputs of this phase will be gvcfs per contig across all samples.

To obtain your list of chromosomes/contigs, first index your .fasta reference genome. Then extract the contig names. Below is an example:

```
samtools reference.fasta # index your fasta first 
cut -f1 reference.fasta.fai > contigs.txt
```

# Phase 4
Phase 4: Run GATK GenotypeGVCFs per contigs. 

Input required for this phase:
   1. GVCF files produced from Phase 3. 
   2. Your reference genome.
   3. A list of the sample names.
   
The outputs of this phase will be genotyped VCFs for each contig that are unfiltered. In the next phase we will apply GATK Best Practice's filters.

# Phase 5
Phase 5: Here, we select the type of variant we want (we are selecting SNPs), but this can be done for say, indels as well. See different options in the [VariantFiltration tools](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration) page. Then, we apply GATK best practices filters to the variants. We use the quality values suggested in the "Hard-filtering" protocol.

Input required for this phase:
   1. VCF files produced from Phase 4. 
   2. Your reference genome.
   
   More information is available on the GATK page:
   
   https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
   
 
   
# Phase 6
Phase 6: This phase is very project-specific and the depth and quality of your dataset will determine the filtering thresholds. Here, we present our script based on our Florida Scrub-Jay data but we encourage researchers to first dig in and visualize your data before deciding on your thresholds. Here is a great guide to get started: https://speciationgenomics.github.io/filtering_vcfs/

Input required for this phase:
   1. GATK-filtered VCFs from phase 5
   2. A text file containing samples that have high % missingness that we'd want to remove FIRST (otherwise these will bias our filtering data-wide, for example if we wanted to filter by mean depth, or mean missingness per site)
   
Here are the filters that we employed as an example:
```
bcftools view -S ^20perc.high.miss.dup.inds.txt passed_allsites_${myContig}.vcf -o allsites_removed_ind.vcf #removed those high missing individuals
bcftools view -i 'QUAL>=30 && FMT/GQ>=20' allsites_removed_ind.vcf -o allsites_GQ20.vcf #sites with minQ & minGQ
bcftools view -e 'MEAN(FMT/DP) < 6 || MEAN(FMT/DP) > 33' allsites_GQ20.vcf -o allsites_DP6_33.vcf #exclude sites depth 6 OR 33, which is half and twice the mean depth of all our specific data
bcftools view -e 'F_MISSING>0.15' allsites_DP6_33.vcf -o allsites_miss_${myContig}.vcf #exclude more than 15% missing sites
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' allsites_miss_${myContig}.vcf > filter1_allsites_${myContig}.vcf ## annotate remaining snp IDs
```


Finally, we will combine all of our per-contigs vcfs that have been filtered into one complete VCF file using ```Picard's GatherVCF``` so that we can run downstream analyses. This requires a list of your per-contig vcfs and their file paths in a .txt to read into the program. 

```
java -Xmx100g -jar /programs/picard-tools-2.8.2/picard.jar GatherVcfs I=GATK_variant_vcfs_2022.list.txt O=/fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2022-July/combined_vcfs/GATK_variants_2022.vcf CREATE_INDEX=true QUIET=true VERBOSITY=WARNING

# Check how many sites you have in total
grep -v "^#" /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2022-July/combined_vcfs/GATK_variants_2022.vcf | wc -l
```







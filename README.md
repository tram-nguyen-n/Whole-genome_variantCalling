# FSJ_variantCalling
Scalable pipeline for variant discovery and quality filtering for paired-end, short-read sequences. This pipeline was developed in collaboration with Dr. Elissa Cosgrove at Cornell University.

Here, we present the pipeline using Florida Scrub-Jay whole-genome sequences. However, this workflow can be adapted to other genomic short-read datasets. This workflow follows the GATK Best Practices protocol and is executed with a SLURM scheduler. Broadly, the pipeline consists of six phases: from alignment to a reference genome, to variant calling using GATK HaplotypeCaller, to additional quality filtering (site depth, missingness, genotyping quality, etc). Note that additional filtering thresholds are specific to our dataset and should be modified accordingly. 

# Dependencies
### Software
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.8
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.36
- [BWA](https://bio-bwa.sourceforge.net/) v0.7.13
- [Samtools](http://www.htslib.org/doc/samtools) v1.17
   - [flagstat](http://www.htslib.org/doc/samtools-flagstat.html)
   - [view](http://www.htslib.org/doc/samtools-view.html)
   - [sort](https://www.htslib.org/doc/samtools-sort.html)


## Phase 1
Phase 1: Run fastqc, trimmomatic, BWA, and sort output bam. 

Input requirements for this phase:
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
  
  4. 

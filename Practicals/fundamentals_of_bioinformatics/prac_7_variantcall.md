# Practical 7 - Variant Calling
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# Introduction
In the last three practicals we performed quality control on Illumina reads and then aligned our reads to the reference genome. We also visualised our alignments with IGV and saw examples of homozygous and herozygous SNPs. 
We are now up to the final step in our analysis - variant calling.  
In the context of the biological question we are investigating, this is the step where we find the genotype of our three Iberian samples at the location of the Eurasian lactase persistence SNP, rs4988235 (also called 13910C>T) and from this, we will be able to predict whether or not these people are lactose tolerant. 


[![Variant calling workflow](https://sbc.shef.ac.uk/wrangling-genomics/img/variant_calling_workflow.png)](https://sbc.shef.ac.uk/wrangling-genomics/04-variant_calling/index.html)

## Overview of today

While sequence alignment is potentially the most important aspect of most high-throughput sequencing pipelines, in whole genome sequencing (WGS) experiments, such as we are simulating here, it is crucial to not only identify where reads have mapped, but regions in which they differ. These regions are called "sequence variants", and can take many forms. Three major types of sequence variation include:

1. Single Nucleotide Polymorphisms or Variants (SNPs/SNVs): single-base pair changes. e.g. A->G
2. Insertion-Deletions (InDels): An insertion or deletion of a region of genomic DNA. e.g. AATA->A
3. Structural Variants: These are large segments of the genome that have been inserted, deleted, rearranged or inverted within the genome.

We are working with short reads and will be focusing on SNPs and InDels. 

## Learning Outcomes

1. How to filter read alignments and why
2. How to deduplicate reads and why
3. How to call and filter variants
4. Know what a VCF is and how to interpret it

## Setup and catchup

Let's activate our software!

```bash
source activate bioinf
```

To catchup, running the code below will get you the required data but 

# **Remove duplicates**

After sorting our read alignments (we did this at the end of the last practical, before IGV visualisation), the next step is marking or removing duplicates. But what are duplicates?

During some sequencing processes (in traditional Illumina sequencing), DNA or cDNA fragments are amplified in a Polymerase Chain Reaction (PCR) using specific adapters and primers.
If the initial unique targets are saturated during this process you can lead to a scenario where replicated fragments are amplified, leading to what we refer to as "Library Duplicates" or "PCR Duplicates".
Additionally, the Illumina machine itself can introduce another duplication effect called an "Optical Duplicate" which is a result of the position of a sequencing cluster on the sequencing flowcell.
Duplicates can be particularly problematic in variant calling because they violate the assumptions of variant identification and have the potential to over-inflate sequencing error.

A few programs exist that can identify these duplicates and mark or remove them from the final BAM file.
To do this, these programs identify sets of read pairs that have the same unclipped alignment start and unclipped alignment end (i.e. the reads start and end at the same alignment position). 
While this works well for paired-end reads, it doesn't work for single-end reads. 
Additionally, count-based sequencing approaches such as ChIP-seq and RNA-seq are generally prone to having high-coverage areas (especially if you have deep sequencing) which may look like duplicates.
Some small non-coding RNAs are also short, so its very likely to have similar alignment starts and ends.
And if you have an experiment that involves a low initial input of DNA or RNA, you're likely to get a high level of PCR duplicates anyway!

If this was an RNA sequencing or ChIP-seq experiment, we would likely run both raw and de-duplicated reads through the same pipeline to compare the results to make sure our duplicates are not effecting our final outcome.
However, because we have paired-end reads and our experiment didn't have a limited amount of DNA input (like we might have if we were working with ancient DNA), so we will simply remove duplicates from all three of our BAM files. 

```bash
# Remove duplicates the samtools way
samtools rmdup 2_align/bam/ERR3241917_sorted.bam 2_align/bam/ERR3241917_dedup.bam

samtools rmdup 2_align/bam/ERR3241921_sorted.bam 2_align/bam/ERR3241921_dedup.bam

samtools rmdup 2_align/bam/ERR3241927_sorted.bam 2_align/bam/ERR3241927_dedup.bam
```

# **Variant Calling**

Variant callers work by counting all reference and alternative alleles at every individual site on the reference genome. 
Because there will be two alleles (e.g. A and B) for each individual reference base (assuming the organism that you are sampling is diploid), then there will be sites which are all reference (AA) or alternate (BB) alleles, which we call a homozygous site. 
If the number of sites is close to 50/50 reference and alternate alleles (AB), we have a heterozygous site. 
We saw examples of homozygous and heterozygous sites in the IGV visualisation and these are the types of features that variant calling algorithms detect. 

There are many variant calling algorithms, which all have advantages and disadvantages in terms of selectivity and sensitivity. 
Many algorithms aim to detect *regions* of the genome where many variants have been called, rather than individual sites, and thus are called *haplotype callers*. 

## What is a VCF?
The standard output of a variant caller is a *variant call format* (VCF) file, a tab-separated file which details information about every sequence variant within the alignment. For information on the VCF file specification see [here](https://samtools.github.io/hts-specs/). The VCF file specification is continuously evolving and is regulated by the [Global Alliance for Genomics & Health](https://www.ga4gh.org/#/fileformats-team).
The VCF file contains everything that you need to know about the sequence variant, the chromosome, position, reference and alternate alleles, the variant quality score and the genotype code (e.g. 0/0, 1/1, 0/1). 
Additionally, the VCF file can be annotated to include information on the region in which a variant was found, such as gene information, whether the variant had an ID (from major databases such as NCBI's dbSNP for example) or whether the variant changed an amino-acid codon or not (synonymous vs non-synonymous sequence variants).

![Variant Call Format](https://github.com/BIG-SA/Intro-NGS-July-2018/raw/master/images/vcfformat.jpg)

## Calling variants

Today we are going to use the haplotype-based caller `freebayes`, which is a [Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment](https://github.com/ekg/freebayes). 

Time to run the variant calling. 
All we need is a reference genome sequence (fasta file), a index of the reference genome (which we did in the last practical), and our BAM files. 

Freebayes assumes a diploid organism by default and doesn't require any pre-filtering of read alignments.  We will be calling variants for all three samples in a single command. 

```bash
# Make sure you're in the project root
cd ~/Practical_alignment

# Run freebayes to create a VCF file. We will call variants for all three samples at once
freebayes \
  -f ref/chr2_sub.fa \
  2_align/bam/ERR3241917_sorted.bam 2_align/bam/ERR3241921_sorted.bam 2_align/bam/ERR3241927_sorted.bam > \
  3_variants/chr2_vars.vcf
```

## Interpreting VCF files

Lets have a look at our VCF file using the command `head`. 
The first part of the file is called the header and it contains all information about the reference sequence, the command that was run, and an explanation of every bit of information that's contained within the *FORMAT* and *INFO* fields of each called variant. 
These lines are denoted by two hash symbols at the beginning of the line ("\#\#"). 
The last line before the start of the variant calls is different to most of the header, as it has one \# and contains the column names for the rest of the file. 

After the header comes the actual variant calls. Use the command below to get just the column names and the first variant in the file. 

```bash
grep "^#CHROM" -A 1 3_variants/chr2_vars.vcf
```

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ERR3241927	ERR3241921	ERR3241917
NC_000002.12	187	.	T	C	581.643	.	AB=0.485714;ABP=3.13438;AC=3;AF=0.5;AN=6;AO=34;CIGAR=1X;DP=70;DPB=70;DPRA=0;EPP=3.26577;EPPR=6.8707;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=3;NUMALT=1;ODDS=27.9738;PAIRED=0.647059;PAIREDR=0.472222;PAO=0;PQA=0;PQR=0;PRO=0;QA=956;QR=1080;RO=36;RPL=19;RPP=4.03217;RPPR=3.0103;RPR=15;RUN=1;SAF=17;SAP=3.0103;SAR=17;SRF=14;SRP=6.8707;SRR=22;TYPE=snp	GT:DP:AD:RO:QR:AO:QA:GL	0/1:16:10,6:10:300:6:180:-11.6811,0,-22.4796	0/1:24:13,11:13:390:11:270:-17.3176,0,-28.1702	0/1:30:13,17:13:390:17:506:-36.8003,0,-26.364

```
### Questions
1. Where is this variant located?
2. What is the reference allele for this variant?
3. What is the alternate allele? 
4. What are the genotypes of the three samples at this location? Are they homozygous or heterozygous?
5. What is the QUAL score?
6. How many variants are there in this file (Hint: each line that's not a header is considered a variant)

<details>
<summary>Answer</summary>
<ul><li>1. On chromosome NC_000002.12 at position 187</li>
<li>2. T</li>
<li>3. C</li>
<li>4. All three are heterozygous (0/1) for this variant.</li>
<li>5. 581.643</li>
<li>6. 22055</li> </ul>
</details>
# **VCF filtering**
Just because a variant is called, does not mean that it is a true positive! We therefore filter variants to keep as many true variants as possible and remove false positives. 

## Filter by QUAL
A very common way to filter variants is by QUAL which is a variant quality score found in the QUAL field. 
QUAL is a phred-scaled quality score for the assertion made in the ALT column and high QUAL scores indicate high confidence calls. 
For example, a QUAL of 30 indicates  a 1 in 1000 chance of this variant not really existing, or 99.9% confidence that the variant does exist.
Variants that pass this filter would be considered "high confidence". 
 

We can use `bcftools filter` to filter for variants that meet certain criteria using an expression filter. 
This means we can either exclude (`-e`) or include (`-i`) variants based on a certain criteria. 

Lets compress our VCF to save space and then exclude all variants that have a QUAL < 30, and then count how many variants meet this criteria. 

```bash
# compress vcf
gzip 3_variants/chr2_vars.vcf

# filter vcf - exclude QUAL<30
bcftools filter -e 'QUAL < 30' 3_variants/chr2_vars.vcf.gz -Oz -o 3_variants/chr2_vars_q30.vcf.gz

# count number of variants that met filtering criteria
zcat 3_variants/chr2_vars_q30.vcf.gz | grep -v "^#" | wc -l

# you could join these two steps together if you don't want to save an intermediate file and just want to count variants passing the filter
bcftools filter -e 'QUAL < 30' 3_variants/chr2_vars.vcf.gz | grep -v "^#" | wc -l
```

- How many variants greater than QUAL 30 do you have? 
- How many variants remain if we exclude variants with QUAL < 20

## Filter by depth - DP

Along with QUAL, variants are also commonly filtered by depth (DP) which can be either total depth of reads across all samples (the INFO DP field) or the depth per sample (the FORMAT DP field). 
Filtering by depth per sample is particularly useful for removing variants called within regions of very high read depth which are often an indication of misalignment or an additional repeat copy in a sample vs the reference. 
Such cases can lead to incorrect calls, which are often extremely confident due to the high depth so won't be removed by QUAL filtering. 
To do this, the value we filter on depends on the average read depth of the sample. If we have an average read depth of around 22x, a good starting point is to remove variants with more than double the average depth, so `FORMAT/DP > 45`. 

Similarly, if you have very low coverage of a variant, it reduces your ability to accurately call a heterozygotic site (especially if you are confident that you sequenced the sample to an adequate depth!). Therefore, we could exclude variants with `INFO/DP<10` (less than 10 reads across all samples aligned to this location). 

Let's re-filter our variants with additional criteria. 

```bash
# filter vcf - exclude q<30 and dp>45
bcftools filter -e 'QUAL<30 | FORMAT/DP>45 | INFO/DP<10' 3_variants/chr2_vars.vcf.gz -Oz -o 3_variants/chr2_vars_q30_dp.vcf.gz

# count number of variants that met filtering criteria
zcat 3_variants/chr2_vars_q30_dp.vcf.gz | grep -v "^#" | wc -l
```

- What do you think the `|` means in the `bcftools filter` command?
- How many variants remained after filtering? 

Now that we have only high-quality variants remaining, we can move on to answering our original question. 
# **Biological Interpretation**

The goal of this analysis was to see how we could use a simple variant calling pipeline to find the genotype of three samples at the site of the Eurasian lactase persistence SNP rs4988235 (also called 13910C>T) and hence infer whether or not these samples were lactose intolerant. 

The Eurasian lactase persistence SNP is located on chromosome 2 at position 135,851,576. In our subset of chromosome 2, this SNP is located at NC_000002.12:2,851,076. 

We can use `bcftools` to filter the VCF file by region or location so let's extract only the Eurasian lactase persistence SNP. 

```bash
# index the filtered vcf so that we can filter by location
bcftools index 3_variants/chr2_vars_q30_dp.vcf.gz

# keep only the variant at NC_000002.12:2851076
bcftools filter -r 'NC_000002.12:2851076' 3_variants/chr2_vars_q30_dp.vcf.gz | less
```

**Questions**
1. What are the reference and alternate alleles? 
2. The MCM6 and LCT genes are located on the reverse strand of the reference genome (see figure below). If the ancestral allele is a C on the reverse strand, what would it be on the forward strand? And similarly, what would the LP allele be on the forward strand if it's a T on the reverse strand? 
3. What are the genotypes of your three samples at this location? 
4. Based on this SNP alone, use the figure below to predict whether or not each of your samples is able to digest lactose as an adult. 

![Lactase Persistence SNP](https://onlinelibrary.wiley.com/cms/asset/98815270-7bb1-4f31-9741-cd3f08d5e1e5/ahg12575-gra-0001.png)

**The genetics of lactase persistence (LP)** from [The molecular basis of lactase persistence: Linking genetics and epigenetics](https://pmc.ncbi.nlm.nih.gov/articles/PMC12336946/)

![[Pasted image 20251024125101.png]]

<details>
<summary>**Answers**</summary>
<ul><li>1. Reference allele is G and alternate is A.</li>
<li>2. Ancestral allele on forward strand is G, LP allele is A.</li>
<li>3. 
ERR3241927: homozygous for alternate allele,
ERR3241921: heterozygous,
ERR3241917: homozygous for reference allele</li>
<li>4. ERR3241927 and ERR3241921 are lactose tolerant while ERR3241917 is lactose intolerant.</li>
</ul>
</details>

# Summary

In this practical, you used a variant calling workflow to genotype the Eurasian lactase persistence SNP in three different samples and used this information to infer lactose tolerance status. 

To do this, you started with raw sequencing data in **FASTQ** format, you first assessed read quality and removed low-quality bases and adapters through **quality control and trimming**. The cleaned reads were then **aligned** to a reference genome, producing **SAM/BAM** files that record where each read maps and how confidently it aligns. After filtering these alignments to remove low-quality and duplicate reads, you performed **variant calling**, generating **VCF** files that list genetic differences between each sample and the reference genome.

While our final interpretation of this analysis focused on a single, well-known variant for teaching purposes, this workflow detected thousands of genetic variants. These variants contain far more information than can be examined manually, including data relevant to health, ancestry, and biological function. Analysing this data is complex and outside the scope of this module but is explored further in another module. 

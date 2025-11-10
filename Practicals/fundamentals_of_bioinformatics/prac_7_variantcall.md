# Practical 7 - Variant Calling
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# Introduction


## Overview of today
In the last three practicals we performed quality control on Illumina reads and then aligned our reads to the reference genome. 
We are now up to the final step in our analysis - variant calling. 

While sequence alignment is potentially the most important aspect of most high-throughput sequencing pipelines, in whole genome sequencing (WGS) experiments, such as we are simulating here, it is crucial to not only identify where reads have mapped, but regions in which they differ. These regions are called "sequence variants", and can take many forms. Three major types of sequence variation include:

1. Single Nucleotide Polymorphisms or Variants (SNPs/SNVs): single-base pair changes. e.g. A->G
2. Insertion-Deletions (InDels): An insertion or deletion of a region of genomic DNA. e.g. AATA->A
3. Structural Variants: These are large segments of the genome that have been inserted, deleted, rearranged or inverted within the genome.

This is where we will finally find the answer to our question! 

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

# Filtering and sorting read alignments

 Section here about choosing params for filtering, and their meanings

Before we start calling variants we will need to sort the alignments. The original file will contain alignments in the order they were found in the original fastq file, so sorting arranges them in genomic order.
```bash
mkdir ~/Practical_7/2_alignedData/sorted_bam
cd ~/Practical_7/2_alignedData
samtools sort bam/SRR2003569_chI.bam -o sorted_bam/SRR2003569_chI.bam
```

This helps the variant caller to run the calling algorithm efficiently and prevents additional time to be allocated to going back and forth along the genome. This is standard for most NGS downstream programs, such as RNA gene quantification and ChIPseq peak calling.

Once we've sorted our alignments, we usually index the file, which allows rapid searching of the file. Running the following command will create an associated bai file which is the index associated with our sorted alignments.

```bash
samtools index sorted_bam/SRR2003569_chI.bam
```


# **Deduplication**

## What are duplicates?
During the sequencing library process, DNA or cDNA fragments are amplified in a Polymerase Chain Reaction (PCR) using specific adapters and primers.
If the initial unique targets are saturated during this process you can lead to a scenario where replicated fragments are amplified, leading to what we refer to as "Library Duplicates" or "PCR Duplicates".
Additionally, the Illumina machine itself can introduce another duplication effect called an "Optical Duplicate" which is a result of the position of a sequencing cluster on the sequencing flowcell.
A few programs exist that can identify these duplicates and remove them from the final BAM file.
To do this, these programs identify sets of reads pairs that have the same unclipped alignment start and unclipped alignment end (i.e. the reads start and end at the same alignment position).

Duplicates can be particularly problematic in variant calling where it violates the assumptions of variant identification and has the potential to over-inflate sequencing error.
However, the catch is that unless you have paired-read sequences (which we do not), its very difficult for you to know that your read is a duplicate in single-end reads because you don't sequence the other end of the fragment.

To add more caveats, count-based sequencing approaches such as ChIP-seq and RNA-seq are generally prone to having high-coverage areas (especially if you have deep sequencing) which may look like duplicates.
Some small non-coding RNAs are also short, so its very likely to have similar alignment starts and ends.
Additionally, if you have an experiment that involves a low initial input of DNA or RNA, you're likely to get a high level of PCR duplicates anyway!

So do we remove duplicates for these approaches?

The answer is:

> ## Evaluate with and without

Generally for RNA sequencing or ChIP-seq experiments, we will run both raw and de-duplicated reads through the same pipeline to compare the results to make sure our duplicates are not effecting our final outcome.

```bash
# Remove duplicates the samtools way
samtools rmdup [SORTED BAM] [SORTED RMDUP BAM]
```

# **Variant Calling**

Variant callers work by counting all reference and alternative alleles at every individual site on the reference genome. 
Because there will be two alleles (e.g. A and B) for each individual reference base (assuming the organism that you are sampling is diploid), then there will be sites which are all reference (AA) or alternate (BB) alleles, which we call a homozygous site. 
If the number of sites is close to 50/50 reference and alternate alleles (AB), we have a heterozygous site.

There are many variant calling algorithms, which all have advantages and disadvantages in terms of selectivity and sensitivity. 
Many algorithms aim to detect *regions* of the genome where many variants have been called, rather than individual sites, and thus are called *haplotype callers*. 

## What is a VCF?
The standard output of a variant caller is a *variant call format* (VCF) file, a tab-separated file which details information about every sequence variant within the alignment. For information on the VCF file specification see [here](https://samtools.github.io/hts-specs/). The VCF file specification is continuously evolving and is regulated by the [Global Alliance for Genomics & Health](https://www.ga4gh.org/#/fileformats-team).
The VCF file contains everything that you need to know about the sequence variant, the chromosome, position, reference and alternate alleles, the variant quality score and the genotype code (e.g. 0/0, 1/1, 0/1). 
Additionally, the VCF file can be annotated to include information on the region in which a variant was found, such as gene information, whether the variant had an ID (from major databases such as NCBI's dbSNP for example) or whether the variant changed an amino-acid codon or not (synonymous vs non-synonymous sequence variants).

![Variant Call Format](https://github.com/BIG-SA/Intro-NGS-July-2018/raw/master/images/vcfformat.jpg)

## Calling variants

Today we are going to use the haplotype-based caller `freebayes`, which is [a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment](https://github.com/ekg/freebayes).

Time to run the variant calling. 
All we need is a reference genome sequence (fasta file), a index of the reference genome (we can do this using `samtools`), and our BAM alignment. 

```bash
# Make sure you're in the project root
cd ~/Practical_alignment
mkdir 3_variants

# Index the reference genome with samtools
samtools faidx ref/chr2_sub.fa

# Run freebayes to create VCF file
freebayes \
  -f ref/chr2_sub.fa \
  2_alignedData/sorted_bam/SRR2003569_chI.bam > \
  2_alignedData/vcf/SRR2003569_chI_1Mb.vcf
```
- What do you think the index is for? 

## Interpreting VCF files

Ok, lets have a look at our VCF file using the command `head`. 
The first part of the file is called the header and it contains all information about the reference sequence, the command that was run, and an explanation of every bit of information that's contained within the *FORMAT* and *INFO* fields of each called variant. 
These lines are denoted by two hash symbols at the beginning of the line ("\#\#"). 
The last line before the start of the variant calls is different to most of the header, as it has one \# and contains the column names for the rest of the file. 
Because we did not specify a name for this sample, the genotype field (the last field of the column line) says "unknown".

```
##fileformat=VCFv4.2
##fileDate=20220904
##source=freeBayes v1.3.2
##reference=genome/chrI.fa
##contig=<ID=I,length=15072434>
##phasing=none
##commandline="freebayes -f genome/chrI.fa --region I:1-1000000 2_alignedData/sorted_bam/SRR2003569_chI.bam"
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
```

After the header comes the actual variant calls, starting from the start of the specified genomic region. 
In our particular file, the header goes for 60 lines, so we may need `sed` to check this section of the file.

```bash
sed -n '61,64p' 2_alignedData/vcf/SRR2003569_chI_1Mb.vcf
```

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  unknown
I       365     .       A       G       5.48143e-08     .       AB=0;ABP=0;AC=0;AF=0;AN=2;AO=2;CIGAR=1X;DP=15;DPB=15;DPRA=0;EPP=7.35324;EPPR=7.18621;GTI=0;LEN=1;MEANALT=1;MQM=4.5;MQMR=11.0769;NS=1;NUMALT=1;ODDS=18.1966;PAIRED=1;PAIREDR=0.692308;PAO=0;PQA=0;PQR=0;PRO=0;QA=81;QR=502;RO=13;RPL=0;RPP=7.35324;RPPR=31.2394;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=9;SRP=7.18621;SRR=4;TYPE=snp GT:DP:AD:RO:QR:AO:QA:GL 0/0:15:13,2:13:502:2:81:0,-3.66059,-12.128
I       371     .       A       G       0.0040121       .       AB=0;ABP=0;AC=0;AF=0;AN=2;AO=3;CIGAR=1X;DP=19;DPB=19;DPRA=0;EPP=9.52472;EPPR=6.62942;GTI=0;LEN=1;MEANALT=2;MQM=20;MQMR=9.13333;NS=1;NUMALT=1;ODDS=6.98659;PAIRED=0.333333;PAIREDR=0.8;PAO=0;PQA=0;PQR=0;PRO=0;QA=112;QR=593;RO=15;RPL=0;RPP=9.52472;RPPR=35.5824;RPR=3;RUN=1;SAF=3;SAP=9.52472;SAR=0;SRF=10;SRP=6.62942;SRR=5;TYPE=snp      GT:DP:AD:RO:QR:AO:QA:GL 0/0:19:15,3:15:593:3:112:0,-0.0395599,-6.91621
I       373     .       G       A       2.09043e-11     .       AB=0;ABP=0;AC=0;AF=0;AN=2;AO=2;CIGAR=1X;DP=20;DPB=20;DPRA=0;EPP=7.35324;EPPR=7.35324;GTI=0;LEN=1;MEANALT=1;MQM=4.5;MQMR=12.3889;NS=1;NUMALT=1;ODDS=26.1012;PAIRED=1;PAIREDR=0.666667;PAO=0;PQA=0;PQR=0;PRO=0;QA=81;QR=656;RO=18;RPL=0;RPP=7.35324;RPPR=42.0968;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=12;SRP=7.35324;SRR=6;TYPE=snp        GT:DP:AD:RO:QR:AO:QA:GL 0/0:20:18,2:18:656:2:81:0,-5.16576,-17.2477

```


### Questions

1. How many variants were called in region ChrI:1-1Mb?
2. The VCF file above shows three variants from the called VCF. What types of variants are they, and are they homozygous or heterozygous?
3. Which INFO field contains information about the variant allele-frequency?



# **VCF filtering**
Just because a variant is called, does not mean that it is a true positive! Each variant called within the file holds a variant quality score (found in the QUAL field). From the VCF format specifications:

QUAL is a phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong). If ALT is ”.” (no variant) then this is -10log_10 p(variant), and if ALT is not ”.” this is -10log_10 p(no variant). High QUAL scores indicate high confidence calls. Although traditionally people use integer phred scores, this field is permitted to be a floating point to enable higher resolution for low confidence calls if desired. (Numeric)

To weed out the low confidence calls in our VCF file we need to filter by QUAL. This can be done using the bcftools program that's included within the samtools suite of tools. All these tools can run on gzip-compressed files which saves a lot of space on your computer.
```bash
gzip SRR2003569_chI_1Mb.vcf
```

Ok lets filter by QUAL. We can do this with the bcftools filter or bcftools view commands which allows you to run an expression filter. This means you can either exclude (-e) or include (-i) variants based on a certain criteria. In our case, lets exclude all variants that have a QUAL < 30.

```bash
bcftools filter -e 'QUAL < 30' SRR2003569_chI_1Mb.vcf.gz -Oz -o SRR2003569_chI_1Mb.sorted.q30.vcf.gz
```

You can pipe to grep and wc to remove the header and count your remaining variants after filtering too

```bash
bcftools filter -e 'QUAL < 30' SRR2003569_chI_1Mb.sorted.q30.vcf.gz | grep -v "^#" | wc -l
```

How many variants greater than QUAL 30 do you have? How about the number of heterozygous variants that have a QUAL>30?

```bash
bcftools filter -i 'QUAL>30 && GT="0/1"' SRR2003569_chI_1Mb.vcf.gz
```
The bcftools view commands gives a lot of additional filtering options.

Questions
- Use the bcftools view or bcftools filter command to count the number of: a. SNPs b. homozygous variants

Depth is also a common filtering characteristic that many people use to remove low confidence variants. If you have low coverage of a variant, it lowers your ability to accurately call a heterozygotic site (especially if you are confident that you sequenced the sample the an adequate depth!). Find the number of SNPs that have a depth that is equal to or greater than 30 and a quality that is greater than 30.


# Biological Interpretation

Discuss likely phenotypes of three individuals. 

# **Visualise long read alignments**





# Read Alignment Practical
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# **Introduction**

In the last two practicals we learnt about quality control. 
We assessed the quality of our raw reads, trimmed them, and then re-assessed their quality to make sure we were happy with the trimming process. 
Now that QC is complete, we move onto the next step: aligning reads to the reference genome.
As a reminder, we will be aligning reads from three Iberian individuals to a small section of human Chromosome 2 to determine the genotype of these individuals at the site of the SNP rs4988235 which predicts lactose tolerance.

[![Variant calling workflow](https://sbc.shef.ac.uk/wrangling-genomics/img/variant_calling_workflow.png)](https://sbc.shef.ac.uk/wrangling-genomics/04-variant_calling/index.html)

## **Overview of today**
A very common approach in bioinformatics is to sequence the genome of an organism (or many organisms) and compare the reads with an existing reference genome.
A reference genome is an accepted representation of an organisms genome that is used by all researchers looking at that particular organism. 
This gives us a shared coordinate system that allows us to document genomic differences in an organised way that is easy to record and communicate. 
In order to make this comparison we must "align" or "map" our reads. 
In read alignment/read mapping, we use software to find the place in the reference that best matches each read (see the figure below). 

![Basic process of aligning reads to a reference genome](https://training.galaxyproject.org/training-material/topics/sequence-analysis/images/mapping/mapping.png)
**Basic process of aligning reads to a reference genome** from [Galaxy mapping training](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html)

This information is documented in a specially formatted text file called a **SAM** file which has the `.sam` extension.
These files tend to be VERY large and so we compress them into **BAM** files which have a `.bam` extension. 
We then filter these alignments in a way that hopefully keeps as many real/true alignments as possible and discards erroneous alignments.
## Learning Outcomes

1. Learn what indexes are for
2. Learn what an alignment file is (SAM and BAM) and understand the contents
3. Run commands to summarise and view mapping statistics from SAM/BAM files
4. How to filter alignments
5. Visualise alignment coverage

## Setup and catchup

Let's activate our `bioinf` conda environment again.

```bash
source activate bioinf
```

We are working in the `~/Practical_alignment` directory again today and will be using the reads that we trimmed in the last practical. 

If you didn't complete the last practical, you may not have all of the required files.
Running the script below will make sure that you are up to date. 

To run the script:
1. Type `nano catchup.sh` to make a new bash script.
2. Paste the entire script below into it and save it by:
	1. Hold down `Ctrl` and type `x` 
	2. Type `y` when you see the message `Save modified buffer?` at the bottom of the screen (this is just asking if you want to save the file)
	3. Press enter when you see `File Name to Write: catchup.sh` at the bottom of the screen. This just means that you are saving the file as `catchup.sh`
3. Run the script with `bash catchup.sh`

```bash
#!/bin/bash

# Sample Variables
SAMPLES=(ERR3241917 ERR3241921 ERR3241927)

# load software
source activate bioinf

# create all directories and move into project directory
mkdir --parents ~/Practical_alignment/{ref,0_raw,1_trim,2_align,3_variants}
mkdir -p ~/Practical_alignment/0_raw/FastQC
mkdir -p ~/Practical_alignment/1_trim/{fastp,FastQC}

cd ~/Practical_alignment

# Get data
## make symlinks
ln -s ~/data/intro_ngs/*.fq.gz 0_raw/
# get reference
cp  ~/data/intro_ngs/chr2_sub.fa ref/

for SAM in "${SAMPLES[@]}";
do
	# Assess raw read quality with fastqc
	fastqc -o 0_raw/FastQC -t 2 0_raw/${SAM}_*.fq.gz

	# Trim with fastp
	fastp --thread 2 \
	-i 0_raw/${SAM}_1.fq.gz \
	-I 0_raw/${SAM}_2.fq.gz \
	-o 1_trim/${SAM}_1.fq.gz \
	-O 1_trim/${SAM}_2.fq.gz \
	--unpaired1 1_trim/${SAM}_1_orphan.fq.gz \
	--unpaired2 1_trim/${SAM}_2_orphan.fq.gz \
	--cut_right \
	--cut_window_size 4 \
	--cut_mean_quality 20 \
	--length_required 75 \
	--html 1_trim/fastp/${SAM}_fastp.html

	# Assess read quality after trimming with fastqc
	fastqc -o 1_trim/FastQC --threads 2 \
	1_trim/${SAM}_{1,2}.fq.gz

done
```


## The reference genome

We will be aligning our reads to a 7Mb portion of the human reference genome GRCh38.p14.  
Let's take a look at the NCBI page for GRCh38.p14 [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) and answer some questions. 
Note that this page includes information on both a RefSeq version and a GenBank version of GRCh38.p14. 
GenBank is a database of all publicly submitted sequences while RefSeq (which stands for 'Reference Sequence') is a subset of GenBank. 
RefSeq provides a single reference genome for each species/strain, often with improved annotations. We are using a portion of the RefSeq reference genome.  

- How long is the human reference genome?
- What is the GC content?
- What is Chromosome 2 called in the RefSeq reference genome?
- How long is Chromosome 2?

The coordinates that describe the portion of GRCh38.p14 that we are analysing are `NC_000002.12:133000000-140000000`.
- NC_000002.12 is the chromosome
- 133000000 is the start position
- 140000000 is the end position

It is called `chr2_sub.fa` and you copied it to your `ref` directory in the first quality control practical. 
Let's take a look. 

```bash
cd ~/Practical_alignment
ls ref/
```

Aligning reads to a reference genome is computationally expensive which means that it takes a long time and a lot of computational resources (CPU and RAM) and this process would be much longer again if we didn't first index the reference genome. 

The first step in aligning reads to a reference genome is to index the reference genome for use by BWA. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment

An index for a reference genome works much the same as an index in a recipe book. 
If you wanted to find all of the recipes that have chocolate in them, you could page through the book and look at every single recipe or you could go to the index, look up chocolate, and get the page numbers for all of the chocolate recipes at once.  
In the case of the reference index, the read alignment tool uses it to look up potential match locations much faster than if it had to start at the beginning of the reference sequence for every read. 

We are using the read alignment tool `bwa` and so we must use `bwa` to index our reference.

```bash
bwa index ref/chr2_sub.fa 
```

This will have created a number of files.

```bash
ls ref/
```

It is not necessary to look through these files as they are specific to the aligner and some of them aren't even human readable. 
## Read alignment
Now that we have trimmed our reads and created an index for our reference, we are ready to start read alignment. 
Let's make a directory for our alignment output files.
```bash
mkdir -p 2_align/{bam,log}
```

The command for read alignment is quite long so we'll talk through it and then run it after. 
This will help you figure out what is going wrong if you get some error messages.

First up here's the command. 

```bash
bwa mem \
  -t 2 \
  ref/chr2_sub.fa \
  1_trim/ERR3241917_1.fq.gz \
  1_trim/ERR3241917_2.fq.gz | \
  samtools view -bhS -F4 - > \
  2_align/bam/ERR3241917_aligned.bam
```

In the first part of the command shown below we use bwa mem to align our trimmed reads to our reference genome. 

```bash
bwa mem \
  -t 2 \
  ref/chr2_sub.fa \
  1_trim/ERR3241917_1.fq.gz \
  1_trim/ERR3241917_2.fq.gz 
```

The default behaviour of bwa mem is to write alignment data to `stdout`, **NOT** to a file. 
You could stream this output to a text file with the ".sam" extension (known as a a SAM file) but SAM files are often very large and so we need a way to compress them.  
To do this, we pipe the alignment data generated by bwa mem to the `samtools` command which converts the plain text SAM format to a binary version, creating a compressed binary SAM file called a BAM file.  

```bash
samtools view -bhS -F4 - 
```

In this context, `samtools` view is the general command that allows the conversion of the SAM to BAM. 
The *globbed* arguments are 1) `-b` [output in binary format]; and 2) `-h` include the file header, followed by the option `-F4` which only include reads with the flag bit `4` set.
(We'll discuss flags in the next section).
The binary output is then written to the file `2_align/bam/ERR3241917_aln.bam` using the `>` symbol.

Here is the command for you to cut and paste:

```bash
bwa mem -t 2 ref/chr2_sub.fa 1_trim/ERR3241917_1.fq.gz 1_trim/ERR3241917_2.fq.gz | samtools view -bhS -F4 - > 2_align/bam/ERR3241917_aln.bam
```
This will take about a minute. 
The output sent to the terminal is normal! 

When it finishes, we can find out information about the alignments using `samtools stats`:

```bash
samtools stats 2_align/bam/ERR3241917_aln.bam > 2_align/log/ERR3241917.stats
```


This is basically the same as the `samtools flagstat` command , but gives additional information.
Take a look at the new file.
```bash
less 2_align/log/ERR3241917.stats
```
### Questions
{:.no_toc}

1. How many reads aligned to our genome?
2. How many reads aligned as a pair?
3. How many aligned as a "proper" pair? ..what is a proper pair??
4. What do you think `inward oriented pairs` and `outward oriented pairs` means?

## SAM (and BAM and CRAM) files


SAM, BAM and CRAM are all different forms of the original SAM format that was defined for holding aligned (or "mapped") high-throughput sequencing data.
SAM stands for Sequence Alignment Map, as described in the [standard specification](http://samtools.github.io/hts-specs/SAMv1.pdf), and was designed to scale to large numbers of reads sequenced at a given time.

The basic structure of each line in the SAM format is depicted in the figure below: 

![](https://us.v-cdn.net/5019796/uploads/editor/f4/uuzmf2cbau1y.png)

SAM files contain a lot of information, with information for every mapped fragment (and sometimes unmapped sequences) being detailed on a single line of text.
Text data generally takes up a large amount of storage space, meaning SAM files are an inefficiant storage format for alignment data.
Instead, storage formats such as BAM and CRAM are often favoured over SAM as they represent the alignment information in a compressed form.
BAM (for Binary Alignment Map) is a lossless compression while CRAM can range from lossless to lossy depending on how much compression you want to achieve (up to very much indeed).
Lossless means that we can completely recover all the data when converting between compression levels, while lossy removes a part of the data that cannot be recovered when converting back to its uncompressed state.
BAMs and CRAMs hold the same information as their SAM equivalent, structured in the same way, however what is different between them is how the files themselves are encoded.

While numbers may vary, generally the file compression that can be achieved from converting the text rich SAM file into the binary BAM version is 1 in 8/10.
Using the default lossless compression in `samtools` (which we will use below), we can almost half the size of a BAM file when converting to a CRAM file.

![File version comparisons](https://www.uppmax.uu.se/digitalAssets/557/c_557912-l_1-k_cram_compression.png)

Most analysis programs that deal with alignments will take SAM and BAM files as input and/or output, and the majority will strictly ask for BAMs as they are more compressed than SAMs.
CRAM files are increasing in popularity and can generally be used with most major programs, with older version containing more limited options for CRAM input.

### Viewing alignments 

To view a SAM, CRAM or BAM file, you can use the [program `samtools`](http://www.htslib.org/).
`samtools` is a very common tool in Bioinformatics and we will be using it frequently in this course.

Lets quickly view our file using the `samtools view` subcommand, which is similar to the command-line tool `cat` in which the file is read to our screen line by line. Make sure you are in your `Project_2` directory.

```bash
samtools view ./data/SRR3096662_Aligned.out.sort.bam
```

As you can probably see, there is a lot of data flashing on your screen.
You can interupt this stream by using/pressing `Ctrl + C` a few times, which should cancel the command.

What you just saw was the alignment information for each read in the `SRR3096662` sample.
Lets use the pipe (`|`) and `head` command to just give us the first five lines of the file so we can start to make sense of the file format.

```text
samtools view ./data/SRR3096662_Aligned.out.sort.bam | head -n5

SRR3096662.22171880 163 1 11680 3 125M = 11681 126 CTGGAGATTCTTATTAGTGATTTGGGCTTGGGCCTGGCCATGTGTATTTTTTTAAATTTCCACTGATGATTTTGCTGCATGGCCGGTGTTGAGAATGACTGCGCAAATTTGCCGGATTTCCTTTG BBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFBF<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFF NH:i:2 HI:i:1 AS:i:244 nM:i:2 MD:Z:28G96 NM:i:1 RG:Z:SRR3096662
SRR3096662.22171880 83 1 11681 3 125M = 11680 -126 TGGAGATTCTTATTAGTGATTTGGGCTTGGGCCTGGCCATGTGTATTTTTTTAAATTTCCACTGATGATTTTGCTGCATGGCCGGTGTTGAGAATGACTGCGCAAATTTGCCGGATTTCCTTTGC FFF<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFBBBB NH:i:2 HI:i:1 AS:i:244 nM:i:2 MD:Z:27G97 NM:i:1 RG:Z:SRR3096662
SRR3096662.15588756 419 1 12009 0 10S115M = 12048 164 GCTTGCTCACGGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCGTTGTTCATCTT BBBBFBFBFFFFFFFFFFFFFFBFFFFFFFFF/FFFFBFFFFFFBB</BF<FFFFFFFFFFFFFFFF/<FFF<FFFFB<FBFB<F/BFFFFBFFBFFFFFBFBFFFFFFFFFFFFFFFFFFFFFF NH:i:5 HI:i:4 AS:i:234 nM:i:2 MD:Z:103A11 NM:i:1 RG:Z:SRR3096662
SRR3096662.15588756 339 1 12048 0 125M = 12009 -164 GCAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCGTTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCATAGGG FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<FFFFFFBFFFFFFFFFFFFFFFFFFFFFFFBFFFFFBFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFBFFBBBB NH:i:5 HI:i:4 AS:i:234 nM:i:2 MD:Z:64A60 NM:i:1 RG:Z:SRR3096662
SRR3096662.17486460 419 1 12174 3 54M385N71M = 12218 549 AAAGATTGGAGGAAAGATGAGTGACAGCATCAACTTCTCTCACAACCTAGGCCAGTGTGTGGTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACT BBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF NH:i:2 HI:i:2 AS:i:241 nM:i:1 MD:Z:24G100 NM:i:1 RG:Z:SRR3096662
```


Each field on each line is followed by a \<TAB> character or what is called a "delimiter".
This just means that the columns in the file are separated by \<TAB> characters, much like a comma-separated file or csv file is delimited by commas.

The specific information of each field is contained below:

| Field | Name | Meaning |
| ---- | ----- | ------- |
| 1 | QNAME | Query template/pair NAME |
| 2 | FLAG | bitwise FLAG (discussed later) |
| 3 | RNAME | Reference sequence (i.e. chromosome) NAME |
| 4 | POS | 1-based leftmost POSition/coordinate of clipped sequence |
| 5 | MAPQ | MAPping Quality (Phred-scaled) |
| 6 | CIGAR | extended CIGAR string |
| 7 | MRNM | Mate Reference sequence NaMe (`=` if same as RNAME) |
| 8 | MPOS | 1-based Mate POSition |
| 9 | TLEN | inferred Template LENgth (insert size) |
| 10 | SEQ | query SEQuence on the same strand as the reference |
| 11 | QUAL | query QUALity (ASCII-33 gives the Phred base quality) |
| 12 | OPT | variable OPTional fields in the format TAG:VTYPE:VALUE |

Notice that each read is considered to be a *query* in the above descriptions, as we a querying the genome to find out where it came from.

Several of these fields contain useful information, so looking the the first few lines, you can see that these reads are mapped in pairs as consecutive entries in the QNAME field are often (but not always) identical.

### Summary statistics


So how do we know that our data is good quality or that our alignment actually worked?
The first basic way of identifying issues is to count how many reads actually mapped!
We can do this easily with the `samtools stats` subcommand, which summarises a lot of quality metrics from the file.
If we extract just the lines starting with "SN" (Summary Numbers), we will get some basic numbers

```bash
samtools stats ./data/SRR3096662_Aligned.out.sort.bam | grep ^SN | cut -f 2-
```

This command will take a while to run, because it summarises all the reference sequences within the CRAM file.
When it finishes, you will see all the summarised information from the file, including aligned reads, how many sequences are found in the header etc...

### SAM headers

While its not actually output in the command that we ran above, an alignment file has a lot of extra information that is commonly printed at the start of a file.
The `samtools view` command actually hides this from you by default, but it contains some very important information such as the names of all the reference genome sequences, whether the file is sorted or not, and the command that was used to create the file in the first place.
To see this info we need to add the `-h` flag.

```bash
samtools view -h ./data/SRR3096662_Aligned.out.sort.bam | less 
```

Page through the output and you should see something like this:

```text
@HD VN:1.4 SO:coordinate
@SQ SN:1 LN:249250621 M5:1b22b98cdeb4a9304cb5d48026a85128 UR:/data/biohub/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta
@SQ SN:2 LN:243199373 M5:a0d9851da00400dec1098a9255ac712e UR:/data/biohub/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta
@SQ SN:3 LN:198022430 M5:fdfd811849cc2fadebc929bb925902e5 UR:/data/biohub/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta
@SQ SN:4 LN:191154276 M5:23dccd106897542ad87d2765d28a19a1 UR:/data/biohub/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta
...
@PG ID:STAR PN:STAR VN:STAR_2.5.2a_modified CL:STAR   --runThreadN 16   --genomeDir /localscratch/Refs/human/hg19_GRCh37d5/star_genome_indices   --readFilesIn ../publicData/SRR3096662_1.fastq.gz   ../publicData/SRR3096662_2.fastq.gz      --readFilesCommand zcat      --outFileNamePrefix ./SRR3096662_   --outSAMtype BAM   SortedByCoordinate      --outSAMattrRGline ID:SRR3096662   LB:library   PL:illumina   PU:machine   SM:hg19      --outFilterType BySJout   --outFilterMismatchNmax 999   --alignSJoverhangMin 8   --alignSJDBoverhangMin 1
@PG ID:samtools PN:samtools PP:STAR VN:1.10 CL:samtools view -@2 -C -T /data/biohub/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta --write-index -o /data/robinson/robertsGroup/180425_Jimmy_eQTL_inputs/crams/SRR3096662_CJM20_Term_Female_Aligned.sortedByCoord.out.cram bams/SRR3096662_CJM20_Term_Female_Aligned.sortedByCoord.out.bam
@PG ID:samtools.1 PN:samtools PP:samtools VN:1.10 CL:samtools view -H SRR3096662_CJM20_Term_Female_Aligned.cram
@RG ID:SRR3096662 LB:library PL:illumina PU:machine SM:hg19
@CO user command line: STAR --genomeDir /localscratch/Refs/human/hg19_GRCh37d5/star_genome_indices --readFilesIn ../publicData/SRR3096662_1.fastq.gz ../publicData/SRR3096662_2.fastq.gz --runThreadN 16 --readFilesCommand zcat --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./SRR3096662_ --outSAMattrRGline ID:SRR3096662 LB:library PL:illumina PU:machine SM:hg19
```

#### Questions

Using the [SAM standard specification](http://samtools.github.io/hts-specs/SAMv1.pdf) and the outputs of the commands shown above, answer the following questions:

1. *Which reference sequence/chromosome did the first read `SRR3096662.22171880` align to and what position did it align to?*
2. *The first two lines of the file contains reads with the same ID (i.e. the first field of the first two lines are `SRR3096662.22171880`). What is a possible reason for this?*
3. *What is the read group ID for the sample?*
4. *What alignment program did we originally use to align the FASTQ data to the reference genome?*
5. *How many mapped reads are contained in our BAM file and what is the average fragment length?*

### CIGAR strings

These give useful information about the type of alignment that has been performed on the read.
In the first few reads we called up earlier, most had the value `..M` where `..` is some number.
These are the perfect Matches, where the sequence has aligned exactly.
The other abbreviations in common use are I (insertion), D (deletion) & S(substitution).


#### Questions
{:.no_toc}

1. What is the interpretation of the first `CIGAR` string in your set of alignments?

### Read Mapping Quality 
Let's run `head` on one of our alignments files again, this time without the header information at the top.

```bash
samtools view ./data/SRR3096662_Aligned.out.sort.bam | head 
```

The 5th field contains the `MAPQ` score which indicates how well the read aligned, and how unique each alignment is.
How this value is calculated can often differ between alignment tools, but primarily the a higher score indicates a better, more unique alignment.

As you have probably learnt from quality scores in sequencing technologies, mapping qualities are projected onto the same quality scale called a "Phred Score".
On sequencing machines, a Phred score is a measure of the probability of error for each individual base pair that was sequenced and therefore is an assigned level of accuracy.
In fact, if you look in the second last field (11th) you will see the ASCII-33 Phred base quality scores.
On the latest sequencing machines that produce FASTQ files (Illumina 1.8+), these scores go from 0 - 41.
The table below defines the error profile of Phred Quality scores as you go higher on the scale.
The higher the score the smaller the probability that you've sequenced an error.

Table: Phred quality scores are logarithmically linked to error probabilities

| Phred Quality Score | Probability of incorrect base call | Base call accuracy |
|---------------------|------------------------------------|--------------------|
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1000 | 99.9% |
| 40 | 1 in 10,000 | 99.99% |
| 50 | 1 in 100,000 | 99.999% |
| 60 | 1 in 1,000,000 | 99.9999% |

Unlike a base quality score, which is the probability that the base is an error, `MAPQ` score is the probability that the sequence read is incorrectly mapped by the aligner.
So by using the `MAPQ` score, we can filter our BAM file to only include the most uniquely, and therefore *most confident*, alignments.
Heng Li (Harvard/MIT), who is the author of a number of alignment algorithms (`bwa` and `minimap2`) and one of the contributing authors for the SAM specification itself, has an [interesting blog post from 2009](http://lh3lh3.users.sourceforge.net/mapuniq.shtml) in which he talks about uniqueness and its relationship to mapping quality:

> Uniqueness was initially introduced to measure the reliability of ungapped short read alignment with a read aligned in full length. It is not a proper concept for generic alignments. For generic alignments, what is much more useful is mapping quality, first introduced in my maq paper. Mapping quality is phred-scaled probability of the alignment being wrong. It unambiguously measures the alignment reliability in a universal way.

There are some issues in how `MAPQ` is implemented in [some alignment algorithms](https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/), but generally it can be used to filter out ambiguously aligned reads that are often found in repetitive regions. (STAR aligner MAPQ values in the range 0-3 use a modified scoring scheme and indicate multi mapping.)

Now that we know what mapping quality is, we can use the `-q` parameter from `samtools view` to filter out any reads that are less than a specific mapping quality.
What mapping quality is a good number to use?
Obviously it is always dependent on what you are trying to do.
If you have very good coverage (i.e. a lot of aligned reads) you can afford to be selective and using a high MAQ such as 30.
For example, in ancient DNA, where DNA fragments are short and therefore likely to have a larger number of repetitive mappings, it is important to filter you data to MAQ => 30.

Lets view the first few alignments that are greater than MAQ30:

```bash
samtools view -q 30 ./data/SRR3096662_Aligned.out.sort.bam | head 
```

As you can see, the first lines have now changed considerably and we only see alignments with >30 quality values.
You can actually go further in filtering mapping quality, and produce a new CRAM file which only contains your high quality alignments.

```bash
samtools view -h -B -q 30 ./data/SRR3096662_Aligned.out.sort.bam -o SRR3096662_Aligned.filtered.bam
```

The additional flags used above are:
`-h`: Print the alignment file with the header
`-B`: Output a BAM file

### SAM Flags
SAM flags are found in the second field of the BAM file and are quite useful pieces of information, however they can be difficult at first look.
Flags can indicate a lot of information that can be used to filter the CRAM file, much like you did with mapping quality.
Head to [http://broadinstitute.github.io/picard/explain-flags.html](https://broadinstitute.github.io/picard/explain-flags.html) to see a helpful description of each flag.

| # | Decimal | Description of read                      |
|---|---------|------------------------------------------|
| 1 | 1       | Read paired                              |
| 2 | 2       | Read mapped in proper pair               |
| 3 | 4       | Read unmapped                            |
| 4 | 8       | Mate unmapped                            |
| 5 | 16      | Read reverse strand                      |
| 6 | 32      | Mate reverse strand                      |
| 7 | 64      | First in pair                            |
| 8 | 128     | Second in pair                           |
| 9 | 256     | Not primary alignment                    |
| 10 | 512    | Read fails platform/vendor quality checks|
| 11 | 1024   | Read is PCR or optical duplicate         |
| 12 | 2048   | Supplementary alignment                  |

Example: for a read with a FLAG value of 163, this is the sum of 128, 32, 2, and 1, which references the 4 descriptions in the table which the read alignment has identified, *"Second in pair - Mate reverse strand - Read mapped in proper pair - Read paired"*.

If we were to identify reads that mapped to the reverse strand, we can use the SAM flag 16.
Then we can use the `-f` parameter to filter our BAM file to only include those reads.

```bash
samtools view -f 16 ./data/SRR3096662_Aligned.out.sort.bam | head
```

You can also do the exact opposite, i.e. identify all reads that are not on the reverse strand, by using the `-F` parameter.

```bash
samtools view -F 16 ./data/SRR3096662_Aligned.out.sort.bam | head
```

Going through a lot of these SAM flags one by one would be fairly tedious, so `samtools` has a subcommand called `flagstat` which counts the number of reads in specific flags.

```bash
samtools flagstat ./data/SRR3096662_Aligned.out.sort.bam
```

### Assess alignment rate and multi-mapping reads

Now that we know what's in a BAM file, how do we assess the quality of the alignment process.
Generally, this is by looking at two metrics:

* how many reads have aligned?
* how many reads multi-mapped?

The primary goal of genome sequence alignment is where you identify the exact position of each read on the reference genome, however it is often the case that a read can map to multiple locations, termed "multi-mapped reads".
In the SAM specifications you can find the full definition of multi-mapping:

> Multiple mapping: The correct placement of a read may be ambiguous, e.g., due to repeats. In this case,
there may be multiple read alignments for the same read. One of these alignments is considered
primary. All the other alignments have the secondary alignment flag set in the SAM records that
represent them. All the SAM records have the same QNAME and the same values for 0x40 and 0x80
flags. Typically the alignment designated primary is the best alignment, but the decision may be
arbitrary.

So each read has a *primary* alignment (i.e. region which the read aligned to which is **usually** the best), and then any subsequent alignment is designated the *secondary* alignment.
Lets take one of the reads from our CRAM file and see whether it is found multiple times

```bash
samtools view  ./data/SRR3096662_Aligned.out.sort.bam  | grep "^SRR3096662.14934677"
```

```
SRR3096662.14934677 99 1 126384 0 125M = 126418 159 TCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGT BBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<FFFFFFF<FBFF NH:i:7 HI:i:1 AS:i:248 nM:i:0 MD:Z:125 NM:i:0 RG:Z:SRR3096662
SRR3096662.14934677 147 1 126418 0 125M = 126384 -159 TTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTA BB<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBB NH:i:7 HI:i:1 AS:i:248 nM:i:0 MD:Z:125 NM:i:0 RG:Z:SRR3096662
SRR3096662.14934677 419 1 336806 0 125M = 336840 159 TAGAGTTCATTTATTTCAACCTGAAGGACTTAACACTTCCTGAATGTCAAATTCAGGGATAAATGGATTTTTTTCAGTTTTAAAAAAAAATCCGGAAATGTCTTAATTTCTCCTTCATTTTTGAA BBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<BB NH:i:7 HI:i:6 AS:i:248 nM:i:0 MD:Z:125 NM:i:0 RG:Z:SRR3096662
SRR3096662.14934677 339 1 336840 0 125M = 336806 -159 ACTTCCTGAATGTCAAATTCAGGGATAAATGGATTTTTTTCAGTTTTAAAAAAAAATCCGGAAATGTCTTAATTTCTCCTTCATTTTTGAAGGATAAGTTTTCCAGCTATATATTTCTCAATTGA FFBF<FFFFFFF<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBB NH:i:7 HI:i:6 AS:i:248 nM:i:0 MD:Z:125 NM:i:0 RG:Z:SRR3096662
SRR3096662.14934677 355 1 652743 0 125M = 652777 159 TCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGT BBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<FFFFFFF<FBFF NH:i:7 HI:i:2 AS:i:248 nM:i:0 MD:Z:125 NM:i:0 RG:Z:SRR3096662
SRR3096662.14934677 403 1 652777 0 125M = 652743 -159 TTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTA BB<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBB NH:i:7 HI:i:2 AS:i:248 nM:i:0 MD:Z:125 NM:i:0 RG:Z:SRR3096662
```

Multi-mapped reads are problematic because we are not confident in their position in the genome, and therefore have a high probability of being erroneous.
This is reflective in the MAPQ score as well.
As you learnt previously, your MAPQ score is the probability that the read is incorrectly mapped, or more importantly, the probability that the read maps *uniquely*.
So if you look at the 5th column of the specific "SRR3096662.14934677" reads, you see that they all have a MAPQ value of zero.
However if you look at the SAM flag field, the reads contain different flags depending on their location.

Additionally, you can check to see how many mappings a read has using the `NH:` tag in the last (12th) field of the line.
That field is optional but many aligners, including the one that we have used i.e. STAR, add the NH tag.
From the STAR manual ():

> The number of loci Nmap a read maps to is given by NH:i: field. Value of 1 corresponds to
unique mappers, while values >1 corresponds to multi-mappers. HI attrbiutes enumerates multiple
alignments of a read starting with 1.

### Questions

Using both `samtools` and other unix tools, identify:

1. *Using the BAM file "SRR3096662_Aligned.out.sort.bam", how many reads differ between `-F 16` and `-f 16`?*
2. *The number of alignments found on Chromosome 1?*
3. *How many aligned reads had a MAPQ greater than 10?*

## Visualise Alignments 


# Prac 7 starts here

## Filtering read alignments

 Section here about choosing params for filtering, and their meanings

## Sorting alignments

Before we start calling variants we will need to **sort the alignments**.
The original file will contain alignments in the order they were found in the original fastq file, so sorting arranges them in *genomic order*.

```bash
mkdir ~/Practical_7/2_alignedData/sorted_bam
cd ~/Practical_7/2_alignedData
samtools sort bam/SRR2003569_chI.bam -o sorted_bam/SRR2003569_chI.bam
```

This helps the variant caller to run the calling algorithm efficiently and prevents additional time to be allocated to going back and forth along the genome. 
This is standard for most NGS downstream programs, such as RNA gene quantification and ChIPseq peak calling.

Once we've sorted our alignments, we usually *index* the file, which allows rapid searching of the file.
Running the following command will create an associated `bai` file which is the index associated with our sorted alignments.

```bash
samtools index sorted_bam/SRR2003569_chI.bam
```

## Deduplication

During the sequencing library process, DNA or cDNA fragments are amplified in a Polymerase Chain Reaction (PCR) using specific adapters and primers.
If the initial unique targets are saturated during this process you can lead to a scenario where replicated fragments are amplified, leading to what we refer to as "Library Duplicates" or "PCR Duplicates".
Additionally, the Illumina machine itself can introduce another duplication effect called an "Optical Duplicate" which is a result of the position of a sequencing cluster on the sequencing flowcell.
A few programs exist that can identify these duplicates and remove them from the final BAM file.
To do this, these programs identify sets of reads pairs that have the same unclipped alignment start and unclipped alignment end (i.e. the reads start and end at the same alignment position).

Duplicates can be particularly problematic in variant calling where it violates the assumptions of variant identification and has the potential to over-inflate sequencing error.
However, the catch is that unless you have paired-read sequences (which we do not), its very difficult for you to know that your read is a duplicate in single-end reads because you don't sequence the other end of the fragment.

To add more caveats, count-based sequencing approaches such as ChIP-seq and RNA-seq are generally prone to having high-coverage areas (especially if you have deep sequencing) which may look like duplcates.
Some small non-coding RNAs are also short, so its very likely to have similar alignment starts and ends.
Additionally, if you have an experiment that involves a low initial input of DNA or RNA, you're likely to get a high level of PCR duplicates anyway!

So do we remove duplicates for these approaches?

The answer is:

> ## Evaluate with and without

Generally for RNA sequencing or ChIP-seq experiments, we will run both raw and de-duplicated reads through the same pipeline to compare the results to make sure our duplicates are not effecting our final outcome.

**For today however, we are going to skip calling duplicates,** but you can use either use `samtools` or `picard` tools to remove them.


## Variant Calling


Variant callers work by counting all reference and alternative alleles at every individual site on the reference genome. 
Because there will be two alleles (e.g. A and B) for each individual reference base (assuming the organism that you are sampling is diploid), then there will be sites which are all reference (AA) or alternate (BB) alleles, which we call a homozygous site. 
If the number of sites is close to 50/50 reference and alternate alleles (AB), we have a heterozygous site.

There are many variant calling algorithms, which all have advantages and disadvantages in terms of selectivity and sensitivity. 
Many algorithms aim to detect *regions* of the genome where many variants have been called, rather than individual sites, and thus are called *haplotype callers*. 
The standard output of a variant caller is a *variant call format* (VCF) file, a tab-separated file which details information about every sequence variant within the alignment. For information on the VCF file specification see [here](https://samtools.github.io/hts-specs/). The VCF file specification is continuously evolving and is regulated by the [Global Alliance for Genomics & Health](https://www.ga4gh.org/#/fileformats-team).
The VCF file contains everything that you need to know about the sequence variant, the chromosome, position, reference and alternate alleles, the variant quality score and the genotype code (e.g. 0/0, 1/1, 0/1). 
Additionally, the VCF file can be annotated to include information on the region in which a variant was found, such as gene information, whether the variant had a ID (from major databases such as NCBI's dbSNP for example) or whether the variant changed an amino-acid codon or not (synonymous vs non-synonymous sequence variants).

![Variant Call Format](https://github.com/BIG-SA/Intro-NGS-July-2018/raw/master/images/vcfformat.jpg)

Today we are going to use the haplotype-based caller `freebayes`, which is [a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment](https://github.com/ekg/freebayes).

Time to run the variant calling. 
All we need is a reference genome sequence (fasta file), a index of the reference genome (we can do this using `samtools`), and our BAM alignment. 
Because variant calling takes a long time to complete, we will only call variants in the first 1Mb of *C. elegans* ChrI to save time. 
To enable us to subset the command, we also need to index the alignment file:

```bash
# Make sure you're in the project root
cd ~/Practical_7
mkdir 2_alignedData/vcf

# Index the reference genome
samtools faidx genome/chrI.fa

# Run freebayes to create VCF file
freebayes \
  -f genome/chrI.fa \
  --region I:1-1000000 \
  2_alignedData/sorted_bam/SRR2003569_chI.bam > \
  2_alignedData/vcf/SRR2003569_chI_1Mb.vcf
```

```bash
freebayes -f genome/chrI.fa --region I:1-1000000 2_alignedData/sorted_bam/SRR2003569_chI.bam > 2_alignedData/vcf/SRR2003569_chI_1Mb.vcf
```

If you're interested in getting all variants from ChrI, run the command without the `--region I:1-1000000` parameter.

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

After the header comes the actual variant calls, starting from the start of the specified genomic region (in our case: ChrI Position 1bp).
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



## VCF filtering
Just because a variant is called, does not mean that it is a true positive! Each variant called within the file holds a variant quality score (found in the QUAL field). From the VCF format specifications:

QUAL is a phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong). If ALT is ”.” (no variant) then this is -10log_10 p(variant), and if ALT is not ”.” this is -10log_10 p(no variant). High QUAL scores indicate high confidence calls. Although traditionally people use integer phred scores, this field is permitted to be a floating point to enable higher resolution for low confidence calls if desired. (Numeric)

To weed out the low confidence calls in our VCF file we need to filter by QUAL. This can be done using the bcftools program that's included within the samtools suite of tools. All these tools can run on gzip-compressed files which saves a lot of space on your computer.

gzip SRR2003569_chI_1Mb.vcf
Ok lets filter by QUAL. We can do this with the bcftools filter or bcftools view commands which allows you to run an expression filter. This means you can either exclude (-e) or include (-i) variants based on a certain criteria. In our case, lets exclude all variants that have a QUAL < 30.

bcftools filter -e 'QUAL < 30' SRR2003569_chI_1Mb.vcf.gz -Oz -o SRR2003569_chI_1Mb.sorted.q30.vcf.gz

# You can pipe to grep and wc to remove the header
#   and count your remaining variants after filtering too
bcftools filter -e 'QUAL < 30' SRR2003569_chI_1Mb.sorted.q30.vcf.gz | grep -v "^#" | wc -l
How many variants greater than QUAL 30 do you have? How about the number of heterozygous variants that have a QUAL>30?

bcftools filter -i 'QUAL>30 && GT="0/1"' SRR2003569_chI_1Mb.vcf.gz
The bcftools view commands gives a lot of additional filtering options.

Questions
Use the bcftools view or bcftools filter command to count the number of: a. SNPs b. homozygous variants

Depth is also a common filtering characteristic that many people use to remove low confidence variants. If you have low coverage of a variant, it lowers your ability to accurately call a heterozygotic site (especially if you are confident that you sequenced the sample the an adequate depth!). Find the number of SNPs that have a depth that is equal to or greater than 30 and a quality that is greater than 30.


## Biological Interpretation

Discuss likely phenotypes of three individuals. 

## Visualise short + long read alignments




# Practical 6 - Read Alignment
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# **Introduction**
This is the third practical in a set of four where we are aligning reads from three Iberian individuals to a small section of human Chromosome 2 and calling variants to determine the genotype of these individuals at the site of the SNP rs4988235. This genotype predicts lactose tolerance/intolerance.

In the last two practicals we learnt about quality control. 
We assessed the quality of our raw reads, trimmed them, and then re-assessed their quality to make sure we were happy with the trimming process. 
Now that QC is complete, we move onto the next step: aligning reads to the reference genome.

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
4. Visualise alignment coverage

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

This 7Mb section of GRCh38.p14 is in a file called `chr2_sub.fa` (as in chromosome 2 subsequence) and ijs locasted in your `ref` directory. 
Let's take a look. 

```bash
cd ~/Practical_alignment
ls ref/
```

- What is the name of this sequence?

The first step in aligning reads to a reference genome is to index the reference genome for use by BWA. 
An index for a reference genome works much the same as an index in a recipe book. 
Indexing allows the alignment tool (`bwa` in this case) to quickly find potential alignment sites for query sequences (reads) in a genome, which saves time during alignment. 

As long as the reference sequence stays the same, we only have to build the index once. However, an index built by `bwa` won't be able to be used by a different read alignment tool. 

```bash
bwa index ref/chr2_sub.fa 
```

This will have created a number of files which are the index.

```bash
ls ref/
```

We don't need to look at the contents of these files as they are specific to the aligner and some of them aren't even human readable. 
## Read alignment
Now that we have trimmed our reads and created an index for our reference, we are ready to start read alignment. 
Let's make a directory for our alignment output files.
```bash
mkdir -p 2_align/{bam,log}
```

Let's look at the `bwa` command for read alignment to see how it works. This will introduce the SAM and BAM file formats. 
SAM files store read mapping/alignment information and BAM files are compressed binary versions of SAM files. 

The command is shown in the figure below with short example filenames. 

![](images/bwa_cmd_explanation.png)

First, BWA mem aligns reads to the reference genome. 
Its default is to write alignment data to `stdout`, not to a file. 
The pipe `|` to captures this stdout (which is in SAM format at this point) to send it to SAMtools instead of just printing it all to the terminal. 
The SAMtools `view` command converts the SAM alignment data to BAM format (a compressed binary version of SAM format) and also throws out any reads with the FLAG 4 set (these are reads that didn't align but we talk about FLAGs a bit later). Finally, the binary output from SAMtools is directed to the new file `aln.bam`.  

Here is the command for you to cut and paste with the appropriate paths:
```bash
bwa mem -t 2 ref/chr2_sub.fa 1_trim/ERR3241917_1.fq.gz 1_trim/ERR3241917_2.fq.gz | samtools view -bh -F4 - > 2_align/bam/ERR3241917_aln.bam
```

This will take about a minute to run and the output sent to the terminal is normal! 

## SAM file format

The output file of the `bwa` command above is a BAM file (a compressed SAM file) which isn't human-readable but is much smaller than a SAM file. 
To see the contents of a BAM file, we use SAMtools to translate the BAM file to SAM format and use `less` to page through it. Remember that `q` will let you exit when you're ready. 

```bash
samtools view -h 2_align/bam/ERR3241917_aln.bam | less
```

SAM stands for Sequence Alignment Map and SAM files describe how individual reads map to a genome. 

The lines beginning with `@` are the header and the rest of the file contains the read alignments. 

In the header, lines beginning with`@SQ` contain reference genome information and lines beginning with `@PG` contain program/software information. 
There are other codes that are often found in SAM headers which you can look up in the [official SAM file documentation](https://samtools.github.io/hts-specs/SAMv1.pdf).  

The column names in the read alignment section are in the table below. 
Each field on each line is followed by a \<TAB> character or what is called a "delimiter".
This just means that the columns in the file are separated by \<TAB> characters, much like a comma-separated file or csv file is delimited by commas.

| Field | Name          | Description                                                                                                     |
| ----- | ------------- | --------------------------------------------------------------------------------------------------------------- |
| 1     | QNAME         | Query/read pair name. This is essentially the first element of the original read identifier                     |
| 2     | FLAG          | bitwise FLAG - see below                                                                                        |
| 3     | RNAME         | Reference sequence chromosome or contig name                                                                    |
| 4     | POS           | 1-based leftmost position/coordinate of clipped sequence in RNAME                                               |
| 5     | MAPQ          | Mapping Quality (Phred-scaled)                                                                                  |
| 6     | CIGAR         | Extended CIGAR string - see below                                                                               |
| 7     | MRNM or MNEXT | Mate Reference sequence name (= if same as QNAME)                                                               |
| 8     | MPOS or PNEXT | 1-based mate position                                                                                           |
| 9     | TLEN          | inferred Template Length (insert size)                                                                          |
| 10    | SEQ           | query sequence on the same strand as the reference (the sequence we aligned). `*` if the sequence is not stored |
| 11    | QUAL          | Query quality (The PHRED scores from the fastq file). If SEQ is `*`, QUAL must also be `*`                      |
| 12    | OPT           | variable optional fields in the format `TAG:TYPE:VALUE`                                                         |
**Note:** "Mate" is referring to a reads pair. ie. The mate of Read1 is Read 2. 

The first read in your alignment file should look something like this: 
❌
```
???
```

Let's talk about the FLAG, MAPQ, and CIGAR columns in more detail. 
### SAM FLAGs

The second column in a SAM file contains what appears to be a normal number but it is actually a _bitwise_ field that contains multiple pieces of information about how a read mapped to the reference. 

The table below contains a set of standardised terms that describe how a read aligned to the reference and gives each term a "Decimal" value. These numbers are determined using a binary system and have the property that all unique combinations of these numbers have a unique value when they are summed.
This means that we can specify any combination of descriptions by simply adding their decimal values to get a single number/integer. [Dave Tang's blog](https://davetang.org/muse/2014/03/06/understanding-bam-flags/) provides a deeper explanation of how this works if you're interested in learning about binary. 

| Decimal | Description of read                       |
| ------- | ----------------------------------------- |
| 1       | Read paired                               |
| 2       | Read mapped in proper pair                |
| 4       | Read unmapped                             |
| 8       | Mate unmapped                             |
| 16      | Read reverse strand                       |
| 32      | Mate reverse strand                       |
| 64      | First in pair                             |
| 128     | Second in pair                            |
| 256     | Not primary alignment                     |
| 512     | Read fails platform/vendor quality checks |
| 1024    | Read is PCR or optical duplicate          |
| 2048    | Supplementary alignment                   |

**Example:** for a read with a FLAG value of 163, this is the sum of 128, 32, 2, and 1, which means that: 
- `128` the read is the second in the pair
- `32 `  the reads mate is on the reverse strand
- `2  ` the read mapped in a proper pair
- `1  ` it is a paired read

You don't have to be able to work this out on your own! 
You can use the [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html) page to find out what  any FLAG value means. 
Head there now to decode the FLAG values 83, 99, 97, and 133. A small sketch might be helpful. 

### MAPQ - mapping quality

The 5th field contains the `MAPQ` score which indicates how well the read aligned, and how unique each alignment is.
How this value is calculated can differ between alignment tools, but primarily, a higher score indicates a better, more unique alignment.
The MAPQ score indicates how likely it is that the read is mapped incorrectly and the likelihood is calculated as below (this is the same as for phred scores). 

*Q =* −10log₁₀*MAPQ*

| Phred Quality Score | Probability of incorrect base call | Base call accuracy |
| ------------------- | ---------------------------------- | ------------------ |
| 10                  | 1 in 10                            | 90%                |
| 20                  | 1 in 100                           | 99%                |
| 30                  | 1 in 1000                          | 99.9%              |
| 40                  | 1 in 10,000                        | 99.99%             |
| 50                  | 1 in 100,000                       | 99.999%            |
| 60                  | 1 in 1,000,000                     | 99.9999%           |

### CIGAR strings
CIGAR strings describe how the individual bases of a read align to the reference using a very simple code of numbers and letters where: 

- `M` - alignment match (can be either a match or mismatch)
- `I` - Insertion to the reference
- `D` - Deletion from the reference
-  = - Sequence match 
- `X`- Sequence mismatch
- `S` - Soft clipping

For example, the CIGAR `8M2I4M1D3M` is interpreted as:
- `8M`  8 alignment matches
- `2I`  2 insertions
- `4M`  4 alignment matches
- `1D`  1 deletion
- `3M`  3 alignment matches

### SAM storage

SAM files contain a lot of information and are often many Gb in size. 
To reduce storage requirements, storage formats such as BAM and CRAM are often favoured over SAM as they represent the alignment information in a compressed form.

BAM (for Binary Alignment Map) is a lossless compression while CRAM can range from lossless to lossy depending on how much compression you want to achieve (up to very much indeed).
Lossless means that we can completely recover all the data when converting between compression levels, while lossy removes a part of the data that cannot be recovered when converting back to its uncompressed state.
BAMs and CRAMs hold the same information as their SAM equivalent, structured in the same way, however what is different between them is how the files themselves are encoded.

Most analysis programs that deal with alignments will take SAM and BAM files as input and/or output, and the majority will strictly ask for BAMs as they are more compressed than SAMs.
CRAM files are increasing in popularity and can generally be used with most major programs, with older versions containing more limited options for CRAM input.

### Summarising alignments

We can get a summary of how our reads aligned using `samtools stats`. It might take a minute to run. 

```bash
samtools stats 2_align/bam/ERR3241917_aln.bam > 2_align/log/ERR3241917.stats
```

You could alternatively use the `samtools flagstat` command which gives slightly less information. 
Have a look at the output using `less`  (`q` to exit).

```bash
less 2_align/log/ERR3241917.stats
```

It tells us that all of the lines beginning with SN are summary information so let's extract just the lines beginning with SN. 
Make sure you quit `less` first with `q`.

```bash
samtools stats 2_align/bam/ERR3241917_aln.bam | grep ^SN | cut -f 2
```

This command will take a while to run, because it summarises all the reference sequences within the BAM file.
When it finishes, you will see all the summarised information from the file, including aligned reads, how many sequences are found in the header etc...

1. How many reads aligned to the genome and what % of reads does it say this is?
2. Is this the same number of reads as were in your trimmed data according to FastQC? What do you think is going on? 
3. How many reads aligned as a pair?
4. How many reads aligned as a "proper" pair? And what is a proper pair?? <- (this might be useful for the assignment)
5. What do you think `inward oriented pairs` and `outward oriented pairs` means?


## Sorting alignments

Before we call variants we will need to **sort the alignments**.
The original file will contain alignments in the order they were found in the original fastq file, so sorting arranges them in *genomic order*.
This helps the variant caller to run the calling algorithm efficiently and prevents additional time to be allocated to going back and forth along the genome. 
This is standard for most NGS downstream programs.

```bash
samtools sort 2_align/bam/ERR3241917_aln.bam -o 2_align/bam/ERR3241917_sorted.bam
```

Once we've sorted our alignments, we usually *index* the file, which allows rapid searching of the file.
Running the following command will create an associated `bai` file which is the index associated with our sorted alignments.

```bash
samtools index 2_align/bam/ERR3241917_sorted.bam
```

We can also delete our unsorted BAM file.  Let's see how big our BAM files are and then delete the unsorted file to reduce storage requirements. 

```bash
ls -lh 2_align/bam/*
rm 2_align/bam/ERR3241917_aln.bam
```

Before we continue, run the read alignment and sorting steps for the two other samples (ERR3241921 and ERR3241927). The code that we used for sample ERR3241917 is shown below. 

```bash
# Align reads
bwa mem -t 2 ref/chr2_sub.fa 1_trim/ERR3241917_1.fq.gz 1_trim/ERR3241917_2.fq.gz | samtools view -bh -F4 - > 2_align/bam/ERR3241917_aln.bam

# Sort
samtools sort 2_align/bam/ERR3241917_aln.bam -o 2_align/bam/ERR3241917_sorted.bam

# index 
samtools index 2_align/bam/ERR3241917_sorted.bam
rm 2_align/bam/ERR3241917_aln.bam
```
## Visualise alignments

We are now going to use IGV to visualise our genome (chr2_sub.fa) and our read alignments. 
We have already sorted and indexed our BAM file (containing read alignments) but need to index the reference sequence as well. 

```bash
samtools faidx ref/chr2_sub.fa
```
Now download the following files to your local computer RStudio's File browser (you can choose a different sample if you want). 
- `~/Practical_alignment/ref/chr2_sub.fa`
- `~/Practical_alignment/ref/chr2_sub.fai`
- `~/Practical_alignment/2_align/bam/ERR3241917_sorted.bam`
- `~/Practical_alignment/2_align/bam/ERR3241917_sorted.bam.bai`

Select 1 file at a time by checking the checkbox and click "More" >> "Export...". Click the "Download" button and save it somewhere obvious.

Visit, [IGV-web](https://igv.org/app/) and load the genome from a `Local File ...` by selecting both the `chr2_sub.fa` and `chr2_sub.fa.fai` files. Once the reference genome is loaded, load a "Track" from a `Local File ...` by selecting both the `.bam` and `.bam.bai` files.

Take some time to work out how to move around, zoom in and out, and then try to answer the questions below:

- The stacked grey arrows represent the reads aligned to the reference genome. What do the coloured vertical bars within the reads indicate?
- Do you see a position in the genome where a coloured vertical bar seems to occur quite frequently? What do you think this is?
- Why does the read coverage drop towards zero at the ends of the reference genome?

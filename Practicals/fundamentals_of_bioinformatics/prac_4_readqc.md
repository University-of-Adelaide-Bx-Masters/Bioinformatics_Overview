# Practical 4 - FASTQ files and Quality Control
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# **Introduction/Background**

## Biological and genetic basis for lactose intolerance 
 
Lactose intolerance affects around 70% of adults worldwide. 
Generally, a healthy newborn baby can digest about 60–70 g of lactose per day (roughly one litre of breast milk) due to the presence of the enzyme lactase in the small intestine. 
However, after weaning, lactase expression typically declines, leading to lactose intolerance in adulthood. 
This is known as lactase non-persistence (LNP) and is the ancestral state for humans.

The lactase gene (LCT) is located on chromosome 2 and is regulated by the MCM6 locus (minichromosome maintenance complex component 6), which lies about 14 kbp upstream. 
Several single nucleotide polymorphisms (SNPs) within the MCM6 locus influence the expression of LCT. 
Certain variants create additional binding sites for transcriptional repressors, reducing the transcription of MCM6 and ultimately lowering lactase production, leading to lactose intolerance. 
Conversely, other SNPs within this region disrupt these repressor sites, allowing lactase expression to persist into adulthood.
This is known as lactase persistence (LP).

With the advent of animal domestication and dairying practices, milk became a reliable and nutrient-rich food source for many populations. 
This shift created strong selective pressure for genetic variants that allowed adults to continue producing lactase and digesting lactose — lactase persistence. 
Because this selection occurred independently in different regions and distinct genetic variants arose in separate populations. 
The distribution of the five most common variants are shown in the figure below. 
Note that all of these variants are located within intron 13 of the MCM6 gene. 

<img src="images/lp_snp_frequencies_by_pop.jpg" alt="SNPs associated with lactose persistence across different populations" height="400">

**SNPs associated with lactose persistence in different populations** from [The molecular basis of lactase persistence: Linking genetics and epigenetics](https://pmc.ncbi.nlm.nih.gov/articles/PMC12336946/)
 
We will be focusing on the Eurasian lactase persistence SNP, [rs4988235](https://asia.ensembl.org/Homo_sapiens/Phenotype/Locations?db=core;name=LACTASE%20PERSISTENCE;ph=3083;r=2:135850576-135851576;v=rs4988235;vdb=variation;vf=89657404), sometimes referred to as 13910C>T as in the figure above.
The A allele enhances activator binding which increases lactase gene expression into adulthood. 

![SNPs conferring lactase persistence](images/lp_snps_in_mcm6.jpg)

**SNPs conferring lactase persistence**
This figure shows a schematic of MCM6 intron 13 lactase persistence enhancer region. 
The light grey represents the genomic sequence in the human reference genome GRCh38 chr2:135,850,966-135,851,196.
The coloured boxes represent transcription factor binding sites and the red lines identify the five SNPs in this region conferring lactase persistence. 
Image from [The molecular basis of lactase persistence: Linking genetics and epigenetics](https://pmc.ncbi.nlm.nih.gov/articles/PMC12336946/)


## Practical Overview

In this practical (and the next three) we will use a simple read alignment and variant calling workflow to determine the genotype of three samples at the site of the rs4988235 SNP.
We will then consider how this relates to the phenotype of lactose tolerance. 
 
The main steps in this workflow are shown in the figure below along with the file types produced by each step. 

[![Variant calling workflow](https://sbc.shef.ac.uk/wrangling-genomics/img/variant_calling_workflow.png)](https://sbc.shef.ac.uk/wrangling-genomics/04-variant_calling/index.html)

The data we will analyse with this workflow is Illumina paired-end reads from three Iberian individuals sequenced as part of the [1000 Genomes project](https://www.coriell.org/1/NHGRI/Collections/1000-Genomes-Project-Collection/1000-Genomes-Project?gad_source=1&gad_campaignid=10942056189&gbraid=0AAAAACRxwMsdRVvA7OauKN189ncoe-14z&gclid=Cj0KCQjwsPzHBhDCARIsALlWNG2QLO7P-lzVqNwqHFEiqk7yXlSRMsX5fLr86aNfAq15Xk-_8Iv5caMaAgmBEALw_wcB)  and the reference genome is a 7Mbp segment (7 million basepairs) from Chromosome 2 in the human reference genome [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). 

The first step in a bioinformatics analysis/workflow is _always_ quality control (QC) and that will be the focus for today. 
This includes checking the quality of raw data, trimming our raw data, and then re-checking quality.
We have included an extra sample for the initial quality assessment so that you see reads of varying quality. 
## Learning Outcomes

1. Gain familiarity with high throughput sequencing data files (FASTQ reads)
2. Learn how to assess the quality of FASTQ reads
3. Learn how to perform adapter and quality trimming
4. Learn simple strategies to reduce memory/storage requirements 

# **Setup**

This practical will again be using RStudio to interact with our VM's. 
See [the first practical](../Bash_Practicals/1_IntroBash.md#rstudio) to remind yourself how to connect. 

All of the code/commands in this practical should be run in the terminal pane. 

## Activate software 

The practicals use an anaconda (`conda`) software environment to provide access to the software you'll need. This is very common practice in bioinformatics. We have set up these environments already so you just need to activate them. 

For today's practical, you will need to activate the `bioinf` conda environment:

```bash
source activate bioinf
```

If this command works properly, it should produce no output. 
**If this command fails and gives you an error message** do the following:

```bash
echo -e "envs_dirs:\n- /apps/conda3/singularity/envs" > ~/.condarc

source activate bioinf
```

Everyones prompt should now have changed to look something like below:

```bash
(bioinf) a1234567@ip-10-255-0-115:/shared/a1234567$
```

The `(bioinf)` prefix lets you know you are in the `bioinf` conda environment, with access to the packages/tools installed in that environment. 
It gives us access to both of the tools (`fastqc` and `fastp`) we need for todays practical. 

## Create directory structure 

Let's create a new directory for todays practical and set up the directory structure for today and the next three practicals. The command `tree` shows the the structure of the `Practical_alignment` directory. 

```bash
mkdir --parents ~/Practical_alignment/{ref,0_raw,1_trim,2_align,3_variants}
# Take a look at the directory structure we created. 
tree Practical_alignment
```

It should look something like below. 

```
Practical_alignment/
├── 0_raw
├── 1_trim
├── 2_align
├── 3_variants
└── ref
```

* *In the `mkdir` command, what did the argument `--parents` do?*
* *In the same command, what was the effect of placing `ref`, `0_raw`, `1_trim` etc.  inside the curly braces?*

## Get data (with symlinks!)

The data we use and create in bioinformatics is often very large and takes up a lot of storage spaces. 
Because data storage is usually limited and is surprisingly expensive, we need to manage our data carefully.
One way we can do that is by using symlinks. 

A symlink (or 'symbolic link') is a shortcut that points to another file or folder. 
It lets you access the original file from a different location without duplicating it and thereby saves disk space.
We will create symlinks to the .fastq files for our analysis and will copy over the reference sequence. 

```bash
cd Practical_alignment

## copy reference genome
cp  ~/data/intro_ngs/chr2_sub.fa ref/

## make symlinks for all .fq.gz files at once
ln -s ~/data/intro_ngs/*.fq.gz 0_raw/

tree .
```
<details>
<summary>Tip! Naming Symlinks</summary>
Symlinks don't have to have the same name as the file they are pointing to but we keep them the same here. 
</details>

The directory structure should now be as below. 

![Directory structure](images/prac4_dir_structure.png)

Notice how the reference sequence is listed just by its name in black while the symlinks are blue with an arrow and the full path to the file they are pointing to is in red.

# **Illumina Sequencing**

In order to analyse our data, we need to understand how it was generated. 
We will be analysing paired-end reads from Illumina, the most commonly used short-read sequencing platform. 
Illumina uses the [Sequencing by Synthesis](https://youtu.be/fCd6B5HRaZ8) method which you will have learned about in the course materials.
More information on Illumina sequencing is available on [their website](https://sapac.illumina.com/systems/sequencing-platforms.html).

<details>
<summary>Sequencing By Synthesis Summary</summary>
<ul>Illumina sequencing uses a sequencing-by-synthesis (SBS) approach, where fluorescently labelled nucleotides are incorporated one base at a time as DNA is copied. First, DNA fragments with special adapters are attached to a flow cell surface and amplified into clusters by bridge amplification, so each cluster represents many copies of the same DNA fragment. During sequencing, DNA polymerase adds one fluorescently labelled nucleotide to each growing strand per cycle. Each base emits a characteristic colour that is imaged by a high-resolution camera, identifying the incorporated base. The fluorescent label and blocking group are then chemically removed so the next base can be added. </ul>

<ul>This process repeats for hundreds of cycles, producing a series of colour images that are computationally converted into a sequence of bases for each cluster. The result is millions to billions of short reads that can then be aligned to a reference genome or assembled de novo for downstream analysis. Sequencing can be performed as single-end (SE), where only one end of each fragment is read, or as paired-end (PE), where both ends are sequenced, providing more information for accurate alignment and detection of structural variation.</ul>
</details>

## Sequencing Template Components

The figure below shows three examples of Illumina templates.  

![](images/illumina_templates.png)

- What happens if the insert fragment is shorter than the read length? 

The figure below shows the different parts of an Illumina sequencing template. 

![](images/seq-template-pe.jpg)

### Adapters
The adapters contain the sequencing primer binding sites (R1 and R2 Primers), index sequences, and the sites that allow library fragments to attach to the flow cell lawn (amplification elements).
There are a limited number of standard Illumina adapter sequences ([detailed here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314)) and so tools are often able to determine which adapter was used automatically.  
### Insert
This is the fragment of DNA that we want to sequence. If we are using barcodes, the barcode is ligated directly to the DNA fragment and is included in the read. 
### Indexes and Barcodes

Indexes and barcodes are similar in that they allow multiple samples to be pooled together in a single sequencing run and later separated by their unique sequence tags. This is called multiplexing. Separating reads by their indexe or barcode into individual samples is called _de-multiplexing_. 
Indexes are not included in either the forward or reverse read and are commonly used in RNA-seq libraries. Indexed samples are generally de-multiplexed by the sequence provider. 

Barcode sequences are included in the read and we have to de-multiplex samples ourselves. Barcodes and indexes can be combined to further multiplex samples onto one sequencing run. 

### Illumina 3' Quality Drop-Off

In general, Illumina sequencing produces highly accurate reads but read quality does tend to diminish towards the 3' end. 
During the bridge-amplification stage, millions of clusters are created on the flowcell. Each cluster comprises of 1,000-2,000 identical copies of the same template. During the process of sequencing, the polymerases attached each of the thousands of copies in a cluster "advance" one base at a time. At the start (5' end) all the polymerases are in perfect sync; generating a bright, clean, consistant light signal for detecting which base was incorporated. However, as time/number of cycles progresses, some polymerases fall behind while some race in front. The polymerases gradually get further and further out of phase with each other. This leads to dimmer and less clear signals for detection, and thus lower quality base-calls.

### Illumina PolyG artifact
Illumina uses only two fluorescent colours in its chemistry to represent the four bases. 
- C  = red
- T  = green
- A = red and green together
- G = no signal

One limitation of this system is that if the signal from a cluster becomes too weak to detect, the instrument interprets the lack of signal as a string of high confidence **G’s**, even if the real bases are different.
This tends to happen more often near the 3' end of reads. 

# **FASTQ file format**

Illumina reads are stored in FASTQ files with the extension `.fq` or `.fastq`. 
These files are plain-text but are often very large so are commonly compressed using `gzip`. The `.gz` extension is added to signify this.
Most modern bioinformatics tools can read `gzip` compressed files you should keep them compressed unless you are using a tool that specifically requires them to be decompresed. 

Let's take a look at the first 4 lines in one of our FASTQ files.

```bash
zcat 0_raw/unknown_R1.fq.gz | head -n 4
```
<details>
<summary>Code Explanation</summary>
<ul><li>'zcat' unzips the file and sends it to stdout (standard output). </li>
<li>The `|` takes standard output from `zcat` and sends it to `head`. </li>
<li>`head -n 4` takes the first 4 lines and sends them to stdout. Because there is nothing after the `head` command to send the output anywhere else, the output is printed to the terminal. </li></ul>
</details>

You should see something like below: 

```
@SRR12313894.2 2 length=76
GNTTCCATGTCGCTGAGTGGAATCTGTTCGTAGGTGGTCGCATAGCGCAGTGCCCTACTATGATCGCTAACTGGAG
+
A!AAFJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJF-AFJ<JJ7AAJJJJJJFJAF7AJJA
```

The four lines are: 
 1. Read identifier, starting with the `@` symbol
 2. Sequence string
 3. A `+` symbol. The read identifier may also immediately follow but this is uncommon.
 4. Quality string

The quality string has a character for each base in the sequence string. Therefore, the length of the sequence string and the quality string should match. Quality values are numbers from `0` to `93` and are often referred to as "Phred" quality scores. To encode the quality scores as a single character, the scores are mapped to the ASCII table:

Standard ASCII Chart - Hex to Decimal code conversion
![Standard ASCII Chart - Hex to Decimal code conversion](https://cdn.shopify.com/s/files/1/1014/5789/files/Standard-ASCII-Table_large.jpg?10669400161723642407) 
From https://www.commfront.com/pages/ascii-chart 

You will see that the first 33 characters (decimal values of 0-32) are all non-printable or white-space (think space, tab, backspace, bell etc). The first printable character is `!` and this has the decimal value of `33`. This character is used to represent a quality value of `0` while `"` has a decimal value of `34` and represents a quality value of `1` (`34-33`). As such these quality scores are said to be Phred+33 encoded and the quality score is simply obtained by substracting 33 from the decimal value of the character in the quality string. Quality score values usually range from 0 (!) to 40 (I) but some files go as high as 41 (J) or 42 (K). 

If you go digging into old Illumina files, you may find quality values which are Phred+64 encoded. That is, a quality value of `0` is represented by `@` which has a decimal value of `64`. However, Phred+33 encoding is the current standard and is often referred to as Illumina 1.9. 

## Phred scores

Phred quality scores give a measure of the confidence the caller has that the sequence base is correct.
To do this, the quality scores are related to the probability of calling an incorrect base through the formula

*Q =* −10log₁₀*P*

where *P* is the probability of calling the incorrect base.
This is more easily seen in the following table:

| Phred Score | Phred+33 Encoded Character | Probability of Incorrect Base Call | Accuracy of Base Call |
| :---------- | :------------------------- | :--------------------------------- | :-------------------- |
| 0           | `!`                        | 1                                  | 0%                    |
| 10          | `+`                        | 10<sup>-1</sup>                    | 90%                   |
| 20          | `5`                        | 10<sup>-2</sup>                    | 99%                   |
| 30          | `?`                        | 10<sup>-3</sup>                    | 99.9%                 |
| 40          | `I`                        | 10<sup>-4</sup>                    | 99.99%                |
|             |                            |                                    |                       |

* In the read we looked at above, the quality string included the characters `A`, `!`, and `J` (as well as some others). Use the ASCII table to determine what Phred score (a number between 0 and 41) these characters represent. Hint: These reads are Phred+33 encoded. 

# **Quality Control**

Quality control generally consists of 3 steps:
1. Check quality of raw data
2. Trim and clean data
3. Check quality of trimmed data

We will be using FastQC to check the quality of our raw and trimmed data and 'fastp' to do the trimming. 
Both of these tools are very commonly used for quality control with Illumina reads.

Let's look at the help file for FastQC to see how it works.

```bash
fastqc -h | less
```

- What does the '-o' option do?
- What does the `-t` option do?

Type `q` to exit when you're finished. 

Let's run FastQC on one of our samples.  

```bash
# create a directory for FastQC output files
mkdir -p ~/Practical_alignment/0_raw/FastQC

# Make sure you're in the right directory
cd ~/Practical_alignment

# Run fastqc on our unknown sample
fastqc -o 0_raw/FastQC/ -t 2 0_raw/unknown_R*.fq.gz
``` 

The above command:

1. Gave both read1 and read2 for a single sample to fastqc using the wildcard `*` in the place of the value 1 or 2.
2. Specified where to write the output (`-o ~/Practical_alignment/0_raw/FastQC`) and
3. Requested two threads (`-t 2`).

<details>
<summary>What are threads?</summary>
<ul>If you had to plant 100 trees and it took you a minute per tree, it would take you 100 minutes to plant all of the trees on your own.
But if you and three friends were doing the job together, it would only take 25 minutes.
These examples are analogous to a single-threaded process and a multi-threaded process in computing.</ul>

<ul>In a multi-threaded process, a big job is split into many smaller jobs that are then run at the same time (in parallel) which gets the job done a lot faster. 
However, some jobs can't be split up in this way because the steps necessary to complete the task are sequential. 
**Your VMs can run two parallel processes (2 threads)** but the University High Performance Computer (HPC) Phoenix can run over 70.
You can tell if a tool can use multiple threads by looking at the documentatiion.</ul>

</details>

FastQC creates an .html report (and a zip file but we don't need that) for each file provided which can be opened in a web browser. 

To view the reports, use the `Files` pane to navigate to `~/Practical_alignment/0_raw/FastQC` and click on one of the html files you find (**If you don't see any html files call a tutor**). Choose `View in Web Browser` and the file will open in your browser.

## Inspecting FastQC reports

Using the Basic Statistics information at the top of the fastqc report, answer the following questions:

- How many reads are in each of the two FASTQ files?
- How many read pairs are there?
- How many reads in total across both files?
- How long are the reads?
- What is the approximate GC content?
- How many Mbp of sequence information is there across both files?
- Coverage (x or fold) is a measure of how many times you would expect each base in the genome to be covered by an aligned base after aligning all the reads. It is the total number of bases in your reads divided by the genome size. If these reads are from a 7Mbp genome, how many fold coverage do we have?

The left hand menu contains a series of clickable links to navigate through the report, with a quick guideline about each section given as a tick, cross or exclamation mark.

Now let's look at some other parts of the report. 

#### Per Base Sequence Quality

Click on the `Per base sequence quality` hyper-link on the left of the page & you will see a boxplot of the QC score distributions for every position in the read.
These are the Phred scores we discussed earlier, and this plot is usually the first one that bioinformaticians will look at for making informed decisions about the overall quality of the data and settings for later stages of the analysis.
The red dashed line indicates the maximum quality score at that position in the read and the blue line indicates the average quality score at that position. 

* *What do you notice about the quality scores as you move from the 5' to the 3' end of the reads?*  Make sure to look at R1 and R2.

#### Per Base Sequence Content

During the preparation of the sequencing library, the genome is randomly fragmented.
As such we would expect to observe relatively flat horizontal lines across this plot.
Demonstrating that we have indeed sampled fragments at random across the genome.
Depending on the GC content of the species, we might see the A and T lines tracking together but separately from the G and C lines.
It is also relatively common to see a drift towards G towards the end of a read.
This is because most modern Illumina sequencing is performed on a "2-colour" system and G is the absence of both colours.
What 
### Overrepresented Sequences

Here we can see any sequence which are more abundant than would be expected

A normal high-throughput library will contain a diverse set of sequences but sometimes we find that some sequences are more common than we would expect. 
This can be for a range of reasons. 
It could mean that a sequence is highly biologically significant, indicate that the library is contaminated, it could be the result of a sequencing artifact. 
Sometimes you will also see sequences here that match the adapters used in the library prep. 
In small RNA libraries where sequences are not randomly fragmented, the same sequence may naturally be present in a significant proportion of the library.

### Adapter Content

For whole genome shotgun sequencing library preps we expect to see little adapter content in the reads.
If there is a significant up-turn in adapter content towards the 3' end, this may indicate the genomic DNA was over-fragmented during library prep.

#### Student Exercise
Run `fastqc` on the sample ERR3241917 and answer the following questions.

```bash
fastqc -o 0_raw/FastQC -t 2 0_raw/ERR3241917_*.fq.gz
```
- What is the maximum quality score shown in the 
- What is the GC content of ERR3241917?
- What might have caused the overrepresented sequence detected in read2 of sample ERR3241917? 
- 


- Basic Statistics
- Per base sequence quality
- Per base sequence content
- Overrepresented sequences

```bash

```
# Adapter and quality trimming

The second step in quality control is adapter removal and quality trimming. 
We need to remove adapters because they don't match any part of the sequenced genome and this will make read alignment more challenging, and low quality stretches of sequence are removed as they are more likely to contain sequencing errors and therefore, may not accurately represent the sequenced genome.
Trimming can also remove some sequencing artifacts, like polyG artifacts.  
Even if the report from FastQC showed that you have good quality data, it is still best-practice to trim.

There are many tools available for trimming.
We are going to be using `fastp` in this practical. The tool `trimmomatic` ([see here](https://github.com/usadellab/Trimmomatic)) would be a good alternative. 

## Fastp
We will use Fastp to trim our reads.  Let's run it and then discuss the code. 
```bash
mkdir -p 1_trim/fastp

fastp --thread 2 -i 0_raw/ERR3241917_1.fq.gz -I 0_raw/ERR3241917_2.fq.gz -o 1_trim/ERR3241917_1.fq.gz -O 1_trim/ERR3241917_2.fq.gz --unpaired1 1_trim/ERR3241917_1_orphans.fq.gz --unpaired2 1_trim/ERR3241917_2_orphans.fq.gz ---cut_right --cut_window_size 4 --cut_mean_quality 20 --length_required 75 --html 1_trim/fastp/ERR3241917_fastp.html
```

The command will take a minute or two to run but it won't tell you that! It looks like it's just hanging/doing nothing. Let's discuss the code while it's running. 
To make the  `fastp` command easier to understand, we have re-formatted it like we might find in a bash script.
```bash
fastp --thread 2 \
-i 0_raw/ERR3241917_1.fq.gz \
-I 0_raw/ERR3241917_2.fq.gz \
-o 1_trim/ERR3241917_1.fq.gz \
-O 1_trim/ERR3241917_2.fq.gz \
--unpaired1 1_trim/ERR3241917_1_orphan.fq.gz \
--unpaired2 1_trim/ERR3241917_2_orphan.fq.gz \
--cut_right \
--cut_window_size 4 \
--cut_mean_quality 20 \
--length_required 75 \
--html 1_trim/fastp/ERR3241917_fastp.html
```

Fastp requires both R1 and R2 fastq files as input and creates four output files. Two files contain the trimmed paired reads (`-o` and `-O` options) and two files contain orphan/unpaired reads. 

Fastp detects and removes adapters from reads by default so we don't have to specifically tell it to do this and the last four parameters/options in the fastp command are for quality trimming. 
Quality trimming is often done using a sliding window. 
This works by moving a small “window” along each read and checking the average quality score of the bases inside that window. In this case, the window is 4 bp wide (`--cut_window_size`) and bases are trimmed if the average quality in the window is less than 20 (`--cut_mean_quality`). 
The `--cut_right` option means that the window starts at the 3' end of the read  and `--length_required 75` indicates that if a read is less than 75bp long after trimming and removing adapters, it should be discarded. 
If a read is discarded but its pair is not, its pair will be sent to the orphan read file. 
## Quality Control post-trimming 

To assess how the trimming has performed, run FastQC across the PAIRED reads output by fastp. 

```bash
mkdir -p ~/Practical_alignment/1_trim/FastQC

fastqc --threads 2 \
1_trim/FastQC/ERR1949188_?.fq.gz
```

Take a look at the various FastQC HTML report files and compare pre and post trimming plots for:

* Per base sequence quality
* Overrepresented sequences
* Adapter Content
# **Building a script**
While we are only analysing three samples in this practical, bioinformaticians often have to analyse many samples at once. 
Instead of manually writing and running a command for each sample, we can build a script to automate it. 

Let's build a script to trim and run fastqc on all three of our human samples.

First, use nano to create a new text file called `run_1_qc.sh` in the `~/Practical_alignment` directory and paste in the script comment skeleton below that shows the main steps in this workflow so far. 
Because the scripts that you write in your assignments will be required to load software, create all of the appropriate directories, obtain data and then run an analysis, we will create this script in the same way.  

First, let's make the directory structure. 

Go back through the practical and find the code you used to do each of these steps (ignore the "Variables" comment for now). 

Add a shebang to the first line of the file and then copy in the `fastp` and `fastqc` commands we used earlier. 
It should look something like below. 
Feel free to add some simple comments. 

```bash
#!/bin/bash

# Variables

# load software

# create directories

# Get data

# Assess raw read quality with fastqc

# Trim with fastp

# Assess read quality after trimming with fastqc

```


```bash
#!/bin/bash

# Variables

# load software
source activate bioinf

# create directories
mkdir --parents ~/Practical_alignment/{ref,0_raw,1_trim,2_align,3_variants}
mkdir -p ~/Practical_alignment/0_raw/FastQC


# Get data
## make symlinks
ln -s ~/data/intro_ngs/*.fq.gz 0_raw/

# Assess raw read quality with fastqc
fastqc -o 0_raw/FastQC -t 2 0_raw/ERR3241917_*.fq.gz

# Trim with fastp

# Assess read quality after trimming with fastqc

```


```bash
#!/bin/bash

## Trim reads
fastp --thread 2 \
-i 0_raw/ERR3241917_1.fq.gz \
-I 0_raw/ERR3241917_2.fq.gz \
-o 1_trim/ERR3241917_1.fq.gz \
-O 1_trim/ERR3241917_2.fq.gz \
--unpaired1 1_trim/ERR3241917_1_orphan.fq.gz \
--unpaired2 1_trim/ERR3241917_2_orphan.fq.gz \
--cut_right \
--cut_window_size 4 \
--cut_mean_quality 20 \
--length_required 75 \
--html 1_trim/fastp/ERR3241917_fastp.html


## Run Fastqc
fastqc --threads 2 \
1_trim/ERR3241917_?.fq.gz
```

If we were to run this script as it is, it would only run on one sample - ERR3241917 - but we want to be able to run it on samples with different names. 
To do that, let's make the sample name into a variable and replace all of the sample names in our code with the variable. Take a look at the script below to see how this might look.  


```bash
#!/bin/bash

# Create a variable called SAMPLE
SAMPLE="ERR3241917"

## Trim reads
fastp --thread 2 \
-i 0_raw/${SAMPLE}_1.fq.gz \
-I 0_raw/${SAMPLE}_2.fq.gz \
-o 1_trim/${SAMPLE}_1.fq.gz \
-O 1_trim/${SAMPLE}_2.fq.gz \
--unpaired1 1_trim/${SAMPLE}_1_orphan.fq.gz \
--unpaired2 1_trim/${SAMPLE}_2_orphan.fq.gz \
--cut_right \
--cut_window_size 4 \
--cut_mean_quality 20 \
--length_required 75 \
--html 1_trim/fastp/${SAMPLE}_fastp.html

## Run Fastqc
fastqc --threads 2 \
1_trim/${SAMPLE}_?.fq.gz
```

For these tools to run, the directories for their output files have to be created first. 

Now change the value of SAMPLE to one you haven't trimmed yet (ERR3241921).

Now save your script and and run it with `bash run_1_qc.sh`. 

```bash
#!/bin/bash

# Create a variable called SAMPLE
SAMPLE="ERR3241917"

## Trim reads
fastp --thread 2 \
-i 0_raw/${SAMPLE}_1.fq.gz \
-I 0_raw/${SAMPLE}_2.fq.gz \
-o 1_trim/${SAMPLE}_1.fq.gz \
-O 1_trim/${SAMPLE}_2.fq.gz \
--unpaired1 1_trim/${SAMPLE}_1_orphan.fq.gz \
--unpaired2 1_trim/${SAMPLE}_2_orphan.fq.gz \
--cut_right \
--cut_window_size 4 \
--cut_mean_quality 20 \
--length_required 75 \
--html 1_trim/fastp/${SAMPLE}_fastp.html

## Run Fastqc
fastqc --threads 2 \
1_trim/${SAMPLE}_?.fq.gz
```

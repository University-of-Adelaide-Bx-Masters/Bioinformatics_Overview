# Introduction to Pathogen Genomics
{:.no_toc}

* TOC
{:toc}

# **1. Introduction**
Imagine the following scenario: 
You are investigating a potential hospital outbreak of *Klebsiella pneumoniae*. As part of this investigation, whole-genome sequencing (WGS) has been performed on 8 *K. pneumoniae* isolates. Your task is to analyse the WGS data to determine which isolates are clonally related.

## 1.1 Practical Overview
In this practical, you will be working with paired-end Illumina sequencing data from 8 bacterial isolates that have been identified as *K. pneumoniae* by microbiological methods. You will perform quality control on the sequencing data and generate a phylogeny to analyse clonal relationships.

## 1.2 Learning Outcomes
1. Gain additional practice in performing read quality control
2. Learn how to perform taxonomic classification using genomic data
3. Learn how to generate and interpret a bacterial core-genome phylogeny


# **2. Setup**

## 2.1 Activate software
For today's practical, you will need to activate the `bioinf` conda environment:

```bash
source activate bioinf
```

## 2.2 Create directory structure
Let's create a new directory for today's practical and create subdirectories that reflect the main steps in our analysis. This will help us stay organised.
```bash
mkdir --parents ~/Practical_pathogens/{ref,db,0_raw,1_trim,2_align,3_phylogeny,4_species}
mkdir -p ~/Practical_pathogens/0_raw/FastQC
mkdir -p ~/Practical_pathogens/1_trim/{fastp,FastQC}
```

## 2.3 Get data
The data for today's practical is located in `~/data/intro_pathogens`. As in previous practicals, we will use symlinks instead of copying large data files.
```bash
cd ~/Practical_pathogens
# create symlinks for all fastq files
ln -s ~/data/intro_pathogens/*.fastq.gz 0_raw/
# create symlink for reference genome
ln -s ~/data/intro_pathogens/CP017934.fasta ref/
# create symlink for kraken database
ln -s ~/data/dbs/kraken/std_8g db/
```

If you run the `tree` command, you can see the structure of all the directories and symlinks you've created. It should look something like this:
```
├── 0_raw
│   ├── FastQC
│   ├── SRR1582858_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582858_30x_1.fastq.gz
│   ├── SRR1582858_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582858_30x_2.fastq.gz
│   ├── SRR1582860_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582860_30x_1.fastq.gz
│   ├── SRR1582860_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582860_30x_2.fastq.gz
│   ├── SRR1582875_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582875_30x_1.fastq.gz
│   ├── SRR1582875_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582875_30x_2.fastq.gz
│   ├── SRR1582878_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582878_30x_1.fastq.gz
│   ├── SRR1582878_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582878_30x_2.fastq.gz
│   ├── SRR1582881_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582881_30x_1.fastq.gz
│   ├── SRR1582881_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582881_30x_2.fastq.gz
│   ├── SRR1582886_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582886_30x_1.fastq.gz
│   ├── SRR1582886_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR1582886_30x_2.fastq.gz
│   ├── SRR2965671_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR2965671_30x_1.fastq.gz
│   ├── SRR2965671_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR2965671_30x_2.fastq.gz
│   ├── SRR2965706_30x_1.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR2965706_30x_1.fastq.gz
│   └── SRR2965706_30x_2.fastq.gz -> /shared//a1078754/data/intro_pathogens/SRR2965706_30x_2.fastq.gz
├── 1_trim
│   ├── fastp
│   └── FastQC
├── 2_align
├── 3_phylogeny
├── 4_species
├── db
│   └── std_8g -> /shared//a1078754/data/dbs/kraken/std_8g
└── ref
    └── CP017934.fasta -> /shared//a1078754/data/intro_pathogens/CP017934.fasta
```


# **3. Read Quality Control**
As in previous practicals, we will use `fastqc` to check the quality of our read data and `fastp` for read trimming.
Since we have 8 samples to process, it will make our code much simpler if we use a `for` loop.
The script below uses steps that you've previously learned about in the Fundamentals of Bioinformatics module. Make sure you understand what each step is doing and ask for help if you're confused.

In order to run the script: 
- Open a file called `run_qc.sh` by typing `nano run_qc.sh`
- Copy the text from the code block below into the `nano` editor
- Save the file and close `nano` by holding down `Ctrl` and typing `x`, then type `y`, then press `Enter`
- To run the script, type `bash run_qc.sh`. It should take around 5 minutes to run. 

```bash
#!/bin/bash

# Load software
source activate bioinf

# List of samples with Illumina data
SAMPLES=(SRR1582858 SRR1582860 SRR1582875 SRR1582878 SRR1582881 SRR1582886 SRR2965671 SRR2965706)

# Loop over each sample
for SAMPLE in "${SAMPLES[@]}";
do
        # Run FastQC on raw data
        fastqc -o 0_raw/FastQC -t 2 0_raw/${SAMPLE}_30x_[12].fastq.gz

        # Perform read trimming, using the same thresholds as in previous practicals
        fastp --thread 2 -i 0_raw/${SAMPLE}_30x_1.fastq.gz -I 0_raw/${SAMPLE}_30x_2.fastq.gz -o 1_trim/${SAMPLE}_1.fastq.gz -O 1_trim/${SAMPLE}_2.fastq.gz --unpaired1 1_trim/${SAMPLE}_1_orphan.fastq.gz --unpaired2 1_trim/${SAMPLE}_2_orphan.fastq.gz --cut_right --cut_window_size 4 --cut_mean_quality 20 --length_required 75 --html 1_trim/${SAMPLE}_fastp.html

        # Run FastQC on trimmed data
        fastqc -o 1_trim/FastQC -t 2 1_trim/${SAMPLE}_[12].fastq.gz
done
```

Now take a look at the FastQC reports. Start with **one sample**  and answer the following questions:
- How many read pairs are present in the raw data and how many are present in the trimmed data?
- What is the total number of bases present in the trimmed reads?
- *K. pneumoniae* has a genome size of approximately 5.5 Mb. Given this, what is the approxiate depth of coverage for this sample?
- Are there any potential quality issues that you think may impact downstream analysis?

Since you have 8 samples, it may take you a lot of time to go through all the FastQC reports in detail. Therefore, instead of carefully inspecting all details, try to skim through the reports to identify any potential quality issues without spending too much time. You will find that most samples have reports that look quite similar to each other so hopefully this shouldn't take you too long.

If you've finished inspecting the FastQC reports and you haven't identified any concerning quality issues, you can move on to the next step. If you're not sure, call over a demonstrator.


# **4. Phylogenetic Analysis**
In order to generate a core genome phylogeny from our WGS data, we need to perform several steps:
- Align our reads to a reference genome
- Perform variant calling
- Generate core SNP alignment
- Build phylogeny

## 4.1 Read Mapping and Variant Calling
You previously learned how to perform read mapping and variant calling in the Fundamentals of Bioinformatics module.
Today we'll use [snippy](https://github.com/tseemann/snippy) to perform these steps. `snippy` uses `bwa mem` for mapping and `Freebayes` for variant calling (like we did in the Fundamentals of Bioinformatics module), however it offers a number of advantages over performing the steps individually. This is because `snippy` is simpler to run, includes important filtering, and produces output in a format that is convenient for generating phylogenies. Because of these advantages, `snippy` is very widely used in pathogen genomics.

Let's start by taking a look at the help message for running snippy:
```bash
snippy --help
```

We can run `snippy` on one of our samples as follows:
```bash
snippy --cpus 2 --ref ref/CP017934.fasta --outdir 2_align/SRR1582858 --R1 1_trim/SRR1582858_1.fastq.gz --R2 1_trim/SRR1582858_2.fastq.gz
```

Take a look at the output files produced:
```bash
ls 2_align/SRR1582858
```

Now you will need to process the remaining samples in the same way by aligning each set of Illumina reads to the same reference. Create a file called `run_snippy.sh` by typing `nano run_snippy.sh` and pasting in the template below. Then fill in the missing steps.
```bash
#!/bin/bash

# Load software

# List of samples with Illumina data

# Loop over each sample

        # Run snippy

```
Once you've finished editing, save the file and close `nano` by holding down `Ctrl` and typing `x`, then type `y`, then press `Enter`.

You can then run your script by typing `bash run_snippy.sh`. This should take around 15 minutes and you will see a lot of output printed to your terminal while it runs. Once this is finished, it's a good idea to check that the expected output files are present.

## 4.2 Generate Core Genome Alignment
Now that we've generated variant calls for each of our samples, the next step is to combine that information to work out which SNPs are "core" and generate an alignment file that can be used as input for tree building. "Core sites" are those positions in the reference genome that are present in *all* the samples. This step will be achieved using `snippy-core`:

```bash
snippy-core --prefix 2_align/core --ref ref/CP017934.fasta 2_align/SRR1582858 2_align/SRR1582860 2_align/SRR1582875 2_align/SRR1582878 2_align/SRR1582881 2_align/SRR1582886 2_align/SRR2965671 2_align/SRR2965706
```

`snippy-core` produces several output files. The one we'll be using in the next step is `core.full.aln`. This is a multi-fasta file that contains the reference genome followed by one sequence for each of our samples. Each of these sequences is the same length as the reference, and the multi-fasta file therefore functions as a multiple sequence alignment. At each position the sequence is: 
- The reference nucleotide if there is no SNP at this position
- The alternative nucleotide if there is a SNP at this position
- '-' if there is a deletion at this position (zero coverage)
- N or n if there is low coverage or a heterozygous or poor quality SNP call


The `core.full.aln` file contains the reference sequence, however we don't want to include this in our phylogeny. Let's generate a modified version that doesn't include the reference.
First, let's find the line numbers for the fasta headers:
```bash
grep -n "^>" 2_align/core.full.aln
```

Your output should look something like this:
```
1:>Reference
89797:>SRR1582858
179593:>SRR1582860
269389:>SRR1582875
359185:>SRR1582878
448981:>SRR1582881
538777:>SRR1582886
628573:>SRR2965671
718369:>SRR2965706
```
This tells us that the reference sequence ends on line 89796, as the header for the second sequence in the file starts on line 89797.
Therefore, we can generate a modified fasta file with no reference sequence by starting at line 89797:
```bash
tail -n +89797 2_align/core.full.aln > 3_phylogeny/core_noref.aln
```

## 4.3 Generate Phylogeny
In order to generate a phylogeny from the core genome alignment, we will use IQ-TREE:

```bash
iqtree -T 2 --mem 16G -s 3_phylogeny/core_noref.aln
```

The output phylogeny is in the file `3_phylogeny/core_noref.aln.treefile`. Let's take a look at the contents of this file:
```bash
cat 3_phylogeny/core_noref.aln.treefile
```

You should see something like this:
```
(SRR1582858:0.0029897898,(SRR1582860:0.0637045832,(SRR1582875:0.0000067158,(SRR1582878:0.0000217506,((SRR1582881:0.0000166650,SRR1582886:0.0000203832):0.0000082761,SRR2965671:0.0000229022):0.0000102206):0.0000238457):0.0023301981):0.0009978756,SRR2965706:0.0031914142);
```

This file is in [Newick format](https://en.wikipedia.org/wiki/Newick_format), which uses parentheses and commas to represent the structure of the tree. This is great for computers, but not easy for humans to understand.
Therefore, we need a visual tool for interpreting our phylogeny.

## 4.4 Visualise Phylogeny
To visualise our phylogeny, we will use the Interactive Tree of Life (iTOL), which is an online tool for tree visualisation.

First, we need to download our tree file. Navigate to the files pane in RStudio and download `core_noref.aln.treefile` to your local machine.

Open [iTOL](https://itol.embl.de/upload.cgi) in a new tab. Select the tree file you have just downloaded and click "Upload". Now your tree looks like a tree!

In order to make the visualisation clearer, we will midpoint root our tree. In the Control panel, select "Advanced", then scroll down to the bottom and select "Midpoint root". This will modify the root of the tree so that it corresponds to the longest branch and is a very simple approach that often makes the structure of a phylogeny much clearer.

If you hover over a branch, the branch length will be displayed. This is measured in substitutions per site. In pathogen genomics, we often refer to the number of substitutions per genome (e.g. if a pair of isolates differ at exactly 5 nucleotide positions within their genomes, we would say they differ by 5 SNPs). To determine the number of SNPs per genome, multiply the branch length by the length of the reference genome (in this case 5.4 Mb).

Questions:
- At the beginning of this practical, you were told that your task was to determine which isolates from this outbreak are clonally related. Can you answer this question now that you have generated the phylogeny?
- You may notice that one isolate is very divergent compared to the other 7 isolates. What is the length of the branch leading to the divergent isolate? Can you think of any possible explanations for this very long branch?


# **5. WGS-based Species Classification**
In pathogen genomics, a very useful quality control step is to perform taxonomic classification of sequencing reads. For example, in the scenario above, we were told that we had 8 *K. pneumoniae* samples. Instead of blindly trusting this information, we should check the sequencing reads for each isolate to make sure that they actually correspond to *K. pneumoniae*.

In this practical, we will use [kraken2](https://ccb.jhu.edu/software/kraken2/) to assign taxonomic labels to sequencing reads using exact k-mer matching, followed by [Bracken](https://ccb.jhu.edu/software/bracken/) to translate the Kraken output into frequency estimates for the constituent taxa.

The `kraken2` database that you will be working with has downsampling of k-mers to give a database size of 8 GB. Without this downsampling, the database would be much larger, and for most applications this would be advantageous, as more k-mers in the database means more accurate classification. However, Kraken requires that the entire database is loaded into memory. As your VMs have only 16 GB memory, an 8 GB database was chosen to enable you to run `kraken2` on your VMs.

## 5.1 Running kraken2

First, let's have a look at the help for `kraken2`:
```bash
kraken2 --help
```

Here is the command for running one of the samples:
```bash
kraken2 --threads 2 --db db/std_8g --output - --report 4_species/SRR1582858.kreport --paired 1_trim/SRR1582858_1.fastq.gz 1_trim/SRR1582858_2.fastq.gz
```

Here is the same command with the options shown on different lines for clarity:
```bash
kraken2 \
    --threads 2 \
    --db db/std_8g \
    --output - \
    --report 4_species/SRR1582858.kreport \
    --paired \
    1_trim/SRR1582858_1.fastq.gz \
    1_trim/SRR1582858_2.fastq.gz
```
This command:
- Uses 2 threads
- Specifies the kraken database as db/std_8g
- Suppresses standard kraken output
- Specifies the file name for the kraken report as 4_species/SRR1582858.kreport
- Specifies that the input reads are paired
- Specifies the file containing the forward reads as 1_trim/SRR1582858_1.fastq.gz
- Specifies the file containing the reverse reads as 1_trim/SRR1582858_2.fastq.gz

The output from running this command is the `kraken2` report, which summarises the number of reads with different taxonomic classifications. We won't look at this directly today, but do feel free to inspect the contents of the file if you're interested.

## 5.2 Running bracken
Bracken (Bayesian Reestimation of Abundance with KrakEN) is a statistical method that estimates the number of reads originating from each species present in a sample. It takes a Kraken report as input and calculates abundance using probabilities derived from the Kraken database that describe how much sequence from each genome is identical to other genomes in the database.

We can run `bracken` as follows:
```bash
bracken -d db/std_8g -r 150 -i 4_species/SRR1582858.kreport -o 4_species/SRR1582858.bracken
```
This command:
- Specifies the kraken2 database that was used as db/std_8g
- Specifies the length of the sequencing reads as 150 bp
- Specifies the input for bracken as the kraken2 report 4_species/SRR1582858.kreport
- Specifies the bracken output location as 4_species/SRR1582858.bracken

Now let's take a look at the `bracken` output:
```bash
cat 4_species/SRR1582858.bracken
```

You should see something like this:
```
name    taxonomy_id     taxonomy_lvl    kraken_assigned_reads   added_reads     new_est_reads   fraction_total_reads
Klebsiella pneumoniae   573     S       37879   341773  379652  0.98927
Klebsiella quasipneumoniae      1463165 S       350     1084    1434    0.00374
Klebsiella variicola    244366  S       168     1368    1536    0.00400
Klebsiella quasivariicola       2026240 S       45      30      75      0.00020
Klebsiella aerogenes    548     S       54      10      64      0.00017
Klebsiella michiganensis        1134687 S       46      71      117     0.00030
Klebsiella africana     2489010 S       33      13      46      0.00012
Klebsiella huaxiensis   2153354 S       26      9       35      0.00009
Klebsiella grimontii    2058152 S       19      86      105     0.00027
Escherichia coli        562     S       153     342     495     0.00129
Salmonella enterica     28901   S       24      181     205     0.00053
```

The important parts to focus on are the first column and the last column. These tell us what the estimated proportion of each species is. In this case, we can see from the second row of the file that 98.927% of reads have been assigned to *K. pneumoniae*. While this is not quite 100%, it is completely normal for a small proportion of reads to be misclassified as closely related species. Therefore, this pattern is exactly what we would expect to see for a pure *K. pneumoniae* sample.

## 5.2 Running all of our samples
Write a script called `run_kraken_bracken.sh` to run the remaining samples through `kraken2` and `bracken`.

Execute your script and then take a look at the `bracken` output files.

Questions:
- Are all of the samples consistent with classification as *K. pneumoniae* isolates?
- Do these results give you any further insight into the divergent sample from section 4?

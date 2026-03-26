# Genome Assembly Quality Control
{:.no_toc}

* TOC
{:toc}

# **1. Introduction**
In previous practicals, you've learned about read quality control. This enables you to identify technical issues that may have occurred during sequencing that could impact your downstream analysis (e.g. poor sequencing quality, contamination with adapter sequences).

However, sometimes data is technically sound, but still problematic. For example, when performing whole-genome sequencing on bacterial isolates:
- Microbiological approaches can misidentify bacterial species
- The isolate may be contaminated
- A sample mixup could mean the isolate that's been sequenced isn't the one that was intended

For these reasons, it is important to check that genomic data matches expectations.
In a previous practical, you learned about using Kraken to perform taxonomic classification on Illumina reads. This is an important step that can identify some types of problematic data.

In this practical, you will learn about additional quality control steps that can be performed on bacterial genome assemblies.

## 1.1 Practical Overview
Scenario: You are a bioinformatician and you've been given Illumina WGS data for 10 *Escherichia coli* and 10 *Streptococcus pneumoniae* isolates. Your task is to perform quality control steps to determine which samples are suitable for downstream analysis.

## 1.2 Learning Outcomes
1. Understand different assembly metrics
2. Learn how assembly metrics can be used to identify problematic data


# **2. Setup**

## 2.1 Activate software
For today's practical, you will need to activate the `bioinf` conda environment:

```bash
source activate bioinf
```

## 2.2 Create directory structure
Let's create a new directory for today's practical and create subdirectories that reflect the main steps in our analysis. This will help us stay organised.
```bash
mkdir -p ~/Practical_assembly_qc/{db,0_raw,1_assembly,2_species}
mkdir -p ~/Practical_assembly_qc/1_assembly/quast
```

## 2.3 Get data
The data for today's practical is located in `~/data/assembly_qc`. As in previous practicals, we will use symlinks instead of copying large data files. Since you will not have time to run 20 genome assemblies during the practical session, we have already generated `spades` assemblies for you. You will need these assemblies as well as the Illumina reads.
```bash
# navigate to working directory
cd ~/Practical_assembly_qc
# create symlinks for fastq files
ln -s ~/data/assembly_qc/*.fastq.gz 0_raw/
# create symlinks for genome assemblies
ln -s ~/data/assembly_qc/*_contigs.fasta 1_assembly
# copy sample lists
cp ~/data/assembly_qc/ecoli_samples.txt .
cp ~/data/assembly_qc/spne_samples.txt .
# create symlink for kraken database
ln -s ~/data/dbs/kraken/std_8g db/
```


# **3. Investigate Assembly Metrics**

## 3.1 *Escherichia coli*
Run `quast` on all *E. coli* assemblies and view the QUAST report. A list of the *E. coli* samples is in the file `ecoli_samples.txt`.

Questions:
- Look at the "# contigs" metric across the 10 samples. Are there any outliers (i.e. do any samples have an unusually high or an unusually low contig number)?
- *E. coli* genomes typically range between 4.5-5.5 Mb. Are any of the samples inconsistent with this expected genome size?
- Expected GC content for an *E. coli* genome is 50-51%. Are any of the samples inconsistent with this range?
- How many samples have at least one suspicious metric identified above?


## 3.2 *Streptococcus pneumoniae*
Run `quast` on all *S. pneumoniae* assemblies and view the QUAST report. A list of the *S. pneumoniae* samples is in the file `spne_samples.txt`.

Questions:
- Are there any samples with an unusually high or unusually low contig number compared to other samples?
- *S. pneumoniae* genomes typically range between 1.9-2.5 Mb. Are any of the samples inconsistent with this expected genome size?
- Expected GC content for a *S. pneumoniae* genome is 39-40%. Are any of the samples inconsistent with this range?
- How many samples have at least one suspicious metric as identified above?
- For the non-suspicious samples, how does contig number and N50 compare to the *E. coli* assemblies?


# **4. Taxonomic classification**
Write and execute a script to run `kraken2` and `bracken` on each of your samples. Investigate each of the bracken output files to determine whether the taxonomic classification is consistent with the expectation of a pure *E. coli* or *S. pneumoniae* isolate.

Questions:
- For the samples with suspicious assembly metrics identified in section 3, does the taxonomic classification explain what's going on?
- Given the taxonomic classification, what do you make of the total assembly length in samples SRR5758738 and SRR6041223?
- Can you explain the assembly metrics for sample SRR7892571?


# **5. Putting it all together**
In this practical, you performed taxonomic classification at the end, after investigating the genome assemblies. This is useful for your learning, because you can see what the impacts are on assemblies when you don't have a pure sample. However, a more typical workflow might look like this:
1. Perform read QC. Exclude QC failures from further analysis.
2. Perform taxonomic classification. Exclude mixed samples and incorrect species from further analysis.
3. Perform genome assembly and analyse assembly metrics. Exclude potentially problematic samples from further analysis (e.g. very high contig number, incorrect genome size, etc).

Another approach that can sometimes be more efficient is to run all steps on all samples, then assess all QC metrics together. This can save time because you can write a single script to do all the analysis steps at once.
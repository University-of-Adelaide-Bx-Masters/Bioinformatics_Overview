
# Common file types in bioinformatics
{:.no_toc}

This document contains brief descriptions and examples of the most commonly encountered file types in bioinformatics. 
There are also links to the official standards/documentations for each file type. 

* TOC
{:toc}

# FASTA

[Official Documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-topics/)  

FASTA format is the simplest format for reporting genomic sequences.
They can store DNA sequences, RNA sequences, or protein sequences of any number and any length.  
If you download a reference genome, it will be in this format (but will likely also be compressed with `gzip`). 
FASTA files (when un-compressed) are human readable. 
### File Extensions

The extensions `.fa`, `.fasta`, and `.fas` are all generic FASTA extensions and they are interchangeable. They can be used for any FASTA file, regardless of its contents. 
There are also a number of other extensions that are used to indicate what type of sequences the file contains. Some of these are shown in the table below. 

| Extension             | Meaning                             |
| --------------------- | ----------------------------------- |
| `.fasta`,`.fas`,`.fa` | Can be used for any FASTA file      |
| `.fna`                | Contains nucleic acid sequence      |
| `.faa`                | Contains amino acid sequence        |
| `.mpfa`               | Contains multiple protein sequences |

### Format

Each sequence in a FASTA file consists of at least two lines.  
1. The first is the header and starts with a `>`
	- Everything up until the first whitespace is the sequence identifier. 
	- Everything after the first white space is the sequence description. Sequence descriptions are not necessary but can be useful.
2. The second is the sequence itself:
	- The sequence can be contained on a single line or take up multiple lines. If the sequence is longer than 80 bp, it is recommended that multiple lines are used. 

### Example

Example of a .fa file containing 3 sequences of different lengths. Sequences longer than 80 bp are wrapped over multiple lines
```fa
>seq1 optional_sequence_description
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTA
>seq2 Sequence description with spaces
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
GAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTT
ATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCT
TTTGGACATTTATTGTCATTCTTACTCCTTTGT
>seq3 [and use brackets]
ATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTG
CATTTGGGAATGTGAGTCTCTTATTGG
```

# FASTQ
[Official Documentation](https://maq.sourceforge.net/fastq.shtml)

FASTQ files store genomic reads and the associated quality scores for each base in the read. When you receive reads from a sequence provider, they will generally be in this form. 
### Extensions
FASTQ files usually have the extension `.fastq`  or `.fq`.
They are commonly compressed with `gzip` and so would have an additional `.gz` extension appended to the end.

### Format
Each individual read in a FASTQ file spans four lines. 
1. The read identifier
2. The sequence of the read itself
3. An alternate line for the identifier but is commonly just a + symbol acting as a placeholder
4. The quality scores for each position along the read as a series of ASCII text characters. 
#### 1. Read Identifier
This line begins with an @ symbol and although there is some variability between different sequencing platforms and software versions, reads from Illumina will traditionally have several components.

An example of a read identifier is shown below and the individual components separated by colons are detailed in the table. 

```
@M02262:117:000000000-AMEE3:1:2105:14127:1629 1:N:0:CGACCTG
```

| Component           | Description                                                                                                         |
| ------------------- | ------------------------------------------------------------------------------------------------------------------- |
|  `@M02262`          | The unique machine ID. This line always begins with an `@` symbol                                                   |
|  `117`              | Run number                                                                                                          |
|  `000000000-AMEE3`  | The flowcell ID                                                                                                     |
|  `1`                | The flowcell lane                                                                                                   |
|  `2105`             | The tile within the flowcell lane                                                                                   |
|  `14127`            | The x-coordinate of the cluster within the tile                                                                     |
|  `1629`             | The y-coordinate of the cluster within the tile                                                                     |
| `1`                 | Indicates that this is the first read in a set of paired-end reads. This would be `2` for the second read in a pair |
| `N`                 | Indicated that this read was NOT marked as a bad quality read by Illumina's own checks                              |
| `0`                 | Value indicating if a read is a control read (0 means it's not)                                                     |
| `CGACCTG`           | The index used when ligating adapters                                                                               |

#### 2. Sequence Read 
This line contains the sequence generated by the sequence platform. 
#### 4. Quality Scores
Quality scores are presented as single ASCII text characters for simple visual alignment with the sequence. 
In the ASCII text system, each character has a numeric value which we can interpret as an integer, and in this context is the quailty score for the corresponding base. Head to the website with a description of these [ASCII Code table](https://en.wikipedia.org/wiki/ASCII#ASCII_printable_code_chart).

The first 31 ASCII characters are non-printable and contain things like end-of-line marks and tab spacings, and note that the first printable character after the space (character 32) is "!" which corresponds to the value 33. 
In short, the values 33-47 are symbols like !, #, $ etc, whereas the values 48-57 are the characters 0-9. Next are some more symbols (including @ for the value 64), with the upper case characters representing the values 65-90 and the lower case letters representing the values 97-122.
### Example
```
@SIM:1:FCX:1:15:6329:1045 1:N:0:2
TCGCACTCAACGCCCTGCATATGACAAGACAGAATC
+
<>;##=><9=AAAAAAAAAA9#:<#<;<<<????#=
```

# GFF 
[Official Documentation]

# GTF
[Official Documentation]

# SAM/BAM
[Official Documentation]

# VCF
[Official Documentation]


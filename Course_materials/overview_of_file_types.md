
# Common file types in bioinformatics
{:.no_toc}

This document contains brief descriptions and examples of commonly encountered file types in bioinformatics. 
There are also links to the official standards/documentation for each file type. 

* TOC
{:toc}

# FASTA

[FASTA Official Documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-topics/)  

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
[FASTQ Official Documentation](https://maq.sourceforge.net/fastq.shtml)

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
#### Read Identifier

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

#### Sequence Read 

This line contains the sequence generated by the sequence platform. 
#### Quality Scores

Quality scores are presented as single ASCII text characters for simple visual alignment with the sequence. 
In the ASCII text system, each character has a numeric value which we can interpret as an integer, and in this context is the quailty score for the corresponding base. Head to the website with a description of these [ASCII Code table](https://en.wikipedia.org/wiki/ASCII#ASCII_printable_code_chart).

The first 31 ASCII characters are non-printable and contain things like end-of-line marks and tab spacings, and note that the first printable character after the space (character 32) is "!" which corresponds to the value 33. 
In short, the values 33-47 are symbols like !, #, $ etc, whereas the values 48-57 are the characters 0-9. Next are some more symbols (including @ for the value 64), with the upper case characters representing the values 65-90 and the lower case letters representing the values 97-122.
### Example
Two Illumina reads are shown below.
```
@HWI-ST999:249:C7M33ACXX:5:1101:3013:2151 1:N:0:TGACCACACTGT
CGCGATAATAATACGCACCGATGACTGGGTGAGAATATTACTTAAGTTCAACAGACTTAAAAATGTTGGGTCCTGGAAAATAATAATCGCCAGC
+
FHHHHHIIJIIJJJJJJJJJJJJJJJJJG@@GIIIGIIJJJJGJJIFH;>EEEDFFFFEEDDD;@@A;;??CDDDDD:>CDDDEECD@BDDB9@
@HWI-ST999:249:C7M33ACXX:5:1101:4346:2125 1:N:0:TGACCACACTGT
GGCACTATCACCGGCGTCTCACGCTTTATGCGCGAACAATCCAAACCGGTGACCATTGTCGGCCTGCAACCGGAAGAGGGCAGCAGCATTCCCG
+
FHHHHHJJJJJJJIJJJJIJJJJIJIIIJIJIJJJGHFDEFEEDEDDDDDDDCDDDDDEEB;BDDDDCDDDDBBBDABBDDBDDDDDBDCDEC<
```

# GFF/GFF3/GTF
[Official Documentation]

# SAM/BAM
[SAM Official Documentation](https://samtools.github.io/hts-specs/SAMv1.pdf) 

SAM stands for Sequence Alignment/Map and SAM files store information that describes how reads map/align to a reference sequence. 
Each line in the file describes the mapping of one read to one location in the reference. 
If a read maps to more than one location, each additional mapping is detailed on a separate line. 
### File Extensions
SAM files will have the `.sam` extension. However, because SAM files tend to be very large, they are often stored as a compressed binary  file - a BAM file - with the extension `.bam`. 

### Format
SAM files begin with a header and header lines start with a `#`. 
The main contents of the file are tab delimited and have at least 9 columns as shown below. 
Columns 10 and 11 may 

| Field | Name  | Description                                                                                     |
| ----- | ----- | ----------------------------------------------------------------------------------------------- |
| 1     | QNAME | Query template/pair NAME. This is essentially the first element of the original read identifier |
| 2     | FLAG  | bitwise FLAG                                                                                    |
| 3     | RNAME | Reference sequence NAME                                                                         |
| 4     | POS   | 1-based leftmost POSition/coordinate of clipped sequence                                        |
| 5     | MAPQ  | mapping Quality (Phred-scaled)                                                                  |
| 6     | CIGAR | extended CIGAR string                                                                           |
| 7     | MRNM  | Mate Reference sequence NaMe ( if same as RNAME)                                                |
| 8     | MPOS  | 1-based Mate POSition                                                                           |
| 9     | TLEN  | inferred Template LENgth (insert size)                                                          |
| 10    | SEQ   | query SEQuence on the same strand as the reference (the sequence we aligned)                    |
| 11    | QUAL  | query QUALity (The PHRED scores from the fastq file)                                            |
| 12    | OPT   | variable optional fields in the format TAG:VTYPE:VALUE                                          |


### Example


# VCF
[VCF Official Documentation](https://samtools.github.io/hts-specs/VCFv4.2.pdf) 

### File Extensions
VCF files use the `.vcf` extension. 

### Format
VCF files begin with a header. All header lines begin with a `#`.
The header includes:

In the main body of the VCF, there are 8 mandatory columns. The  

| Column | Name    | Description                                                                                                                                                                                                      |
| ------ | ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1      | CHROM   | Name of the sequence on which the variation is located (usually a chromosome)                                                                                                                                    |
| 2      | POS     | The 1-based position of the variation on the given sequence.                                                                                                                                                     |
| 3      | ID      | The identifier of the variation, e.g. a [dbSNP](https://en.wikipedia.org/wiki/DbSNP "DbSNP") rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.    |
| 4      | REF     | The reference base (or bases in the case of an [indel](https://en.wikipedia.org/wiki/Indel "Indel")) at the given position on the given reference sequence.                                                      |
| 5      | ALT     | The list of alternative [alleles](https://en.wikipedia.org/wiki/Alleles "Alleles") at this position.                                                                                                             |
| 6      | QUAL    | A quality score associated with the inference of the given alleles.                                                                                                                                              |
| 7      | FILTER  | A flag indicating which of a given set of filters the variation has failed or PASS if all the filters were passed successfully.                                                                                  |
| 8      | INFO    | An extensible list of key-value pairs (fields) describing the variation. See below for some common fields. Multiple fields are separated by semicolons with optional values in the format: `<key>=<data>[,data]` |
| 9      | FORMAT  | An (optional) extensible list of fields for describing the samples. See below for some common fields.                                                                                                            |
| +      | SAMPLEs | For each (optional) sample described in the file, values are given for the fields listed in FORMAT                                                                                                               |
 
### Example
```
##fileformat=VCFv4.2 
##fileDate=20090805 
##source=myImputationProgramV3.1 ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta 
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial 


#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003 
20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,. 20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3 
20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4 20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2 
20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3
```
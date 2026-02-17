
# Common file types in bioinformatics
{:.no_toc}

This document contains brief descriptions and examples of commonly encountered file types in bioinformatics. 
There are also links to the official standards/documentation for each file type. 

All files described here are plain-text files when not compressed.

* TOC
{:toc}

# FASTA

- [FASTA Official Documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-topics/)  

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
- [FASTQ Official Documentation](https://maq.sourceforge.net/fastq.shtml)

FASTQ files store genomic reads and the associated quality scores for each base in the read. When you receive reads from a sequence provider, they will generally be in this form. 

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

# GFF/GTF

- [GFF3 Documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) 

GFF stands for General Feature Format. 
The GFF filetype has been around in some form since 1997 but has morphed into a multitude of different flavours since then (see the [AGAT documentation](https://agat.readthedocs.io/en/latest/gxf.html) for some history). 

The most common are GTF (also called GFF2) and GFF3.  
These file types are all tab delimited with at least 9 columns and describe genomic features. 

GTF files were developed from the original GFF format to better describe complex gene structures, and GFF3 was developed in an effort to standardise the GFF class in a way that was both backwards compatible and flexible. 

Extensions you might encounter include `.gff`, `.gtf`, `.gff2`, and `.gff3`. 
### Format
The version of your GFF file (or GTF) will hopefully be located in the first line of the file. 

All GFF-type files are tab delimited and contain 9 fields. 
Most of these fields are the same, regardless of which type of file you have. 
Below are the column names according to the GFF3 specification. 

| Column | Name      | Description                                                                                                                                                                                                                                                                                                 |
| ------ | --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1      | seqid     | The name of the sequence where the feature is located.                                                                                                                                                                                                                                                      |
| 2      | source    | The algorithm or procedure that generated the feature. This is typically the name of a software or database.                                                                                                                                                                                                |
| 3      | type      | The feature type. For example, `gene` or `exon`, etc.                                                                                                                                                                                                                                                       |
| 4      | start     | Genomic start of the feature with a 1-base offset.                                                                                                                                                                                                                                                          |
| 5      | end       | Genomic end of the feature with a 1-base offset.                                                                                                                                                                                                                                                            |
| 6      | score     | Numeric value that generally indicates the confidence of the source in the annotated feature. A value of "." (a dot) is used to define a null value.                                                                                                                                                        |
| 7      | strand    | Single character that indicates the [strand](https://en.wikipedia.org/wiki/Sense_\(molecular_biology\)#DNA_sense "Sense (molecular biology)") of the feature. This can be "+" (positive, or 5'->3'), "-", (negative, or 3'->5'), "." (undetermined), or "?" for features with relevant but unknown strands. |
| 8      | phase     | Phase of CDS features; it can be either one of 0, 1, 2 (for CDS features) or "." (for everything else).                                                                                                                                                                                                     |
| 9      | attribute | A list of feature atributes in the format `tag=value`. Multiple `tag=value` pairs are separated by semicolons. Some tags with pre-defined meanings are defined in the GFF3 documentation.                                                                                                                   |

The main differences between GFF and GTF files are in column 3 and column 9. 
In column 3, GTF files are constrained to less than 10 possible feature types while GFF files are constrained to 2278 possible feature types. 
In column 9, the structure differs slightly. 
### Example

 A `.gff3` file describing a gene called 'EDEN' with multiple transcripts. Note how the 9th column contains an ID for each feature and also details its parent feature (ie. the parent feature of mRNA00001 is gene00001). This is how individual features are linked together.  

```gff3
##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
ctg123 . TF_binding_site 1000  1012  .  +  .  ID=tfbs00001;Parent=gene00001
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123 . mRNA            1300  9000  .  +  .  ID=mRNA00003;Parent=gene00001;Name=EDEN.3
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00003
ctg123 . exon            1050  1500  .  +  .  ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123 . exon            3000  3902  .  +  .  ID=exon00003;Parent=mRNA00001,mRNA00003
ctg123 . exon            5000  5500  .  +  .  ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123 . exon            7000  9000  .  +  .  ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123 . CDS             1201  1500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             3000  3902  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             5000  5500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             7000  7600  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             1201  1500  .  +  0  ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123 . CDS             5000  5500  .  +  0  ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123 . CDS             7000  7600  .  +  0  ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123 . CDS             3301  3902  .  +  0  ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123 . CDS             5000  5500  .  +  1  ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123 . CDS             7000  7600  .  +  1  ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123 . CDS             3391  3902  .  +  0  ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
ctg123 . CDS             5000  5500  .  +  1  ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
ctg123 . CDS             7000  7600  .  +  1  ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
```

# SAM/BAM

- [SAM Official Documentation](https://samtools.github.io/hts-specs/SAMv1.pdf) 
- [SAM Flag decoder](https://broadinstitute.github.io/picard/explain-flags.html) 

SAM stands for Sequence Alignment/Map and SAM files store information that describes how reads map/align to a reference sequence. 
Each line in the file describes the mapping of one read to one location in the reference. 
If a read maps to more than one location, each additional mapping is detailed on a separate line. 

SAM files will have the `.sam` extension. 
SAM files are very large and are usually stored as a compressed binary file - a BAM file - with the extension `.bam`. 
They may also be further compressed as a CRAM with the extension `.cram`. 
### Format

SAM files begin with a header and all header lines start with `@`. 
The main contents of the file are tab delimited and have at least 11 columns as shown below. 

| Field | Name          | Description                                                                                                     |
| ----- | ------------- | --------------------------------------------------------------------------------------------------------------- |
| 1     | QNAME         | Query template/pair NAME. This is essentially the first element of the original read identifier                 |
| 2     | FLAG          | bitwise FLAG                                                                                                    |
| 3     | RNAME         | Reference sequence chromosome or contig NAME                                                                    |
| 4     | POS           | 1-based leftmost POSition/coordinate of clipped sequence                                                        |
| 5     | MAPQ          | MAPping Quality (Phred-scaled)                                                                                  |
| 6     | CIGAR         | extended CIGAR string                                                                                           |
| 7     | MRNM or MNEXT | Mate Reference sequence NaMe ( if same as RNAME)                                                                |
| 8     | MPOS or PNEXT | 1-based Mate POSition                                                                                           |
| 9     | TLEN          | inferred Template LENgth (insert size)                                                                          |
| 10    | SEQ           | query SEQuence on the same strand as the reference (the sequence we aligned). A `*` if the sequence  not stored |
| 11    | QUAL          | query QUALity (The PHRED scores from the fastq file). If SEQ is `*`, QUAL must also be `*`                      |
| 12    | OPT           | variable optional fields in the format `TAG:TYPE:VALUE`                                                         |

#### Field 2 - SAM Flags

Summary of SAM flags. 

| #   | Decimal | Description of read                       |
| --- | ------- | ----------------------------------------- |
| 1   | 1       | Read paired                               |
| 2   | 2       | Read mapped in proper pair                |
| 3   | 4       | Read unmapped                             |
| 4   | 8       | Mate unmapped                             |
| 5   | 16      | Read reverse strand                       |
| 6   | 32      | Mate reverse strand                       |
| 7   | 64      | First in pair                             |
| 8   | 128     | Second in pair                            |
| 9   | 256     | Not primary alignment                     |
| 10  | 512     | Read fails platform/vendor quality checks |
| 11  | 1024    | Read is PCR or optical duplicate          |
| 12  | 2048    | Supplementary alignment                   |

Example: 
Let's interpret the FLAG value 163. 

This is the sum of 1, 2, 32, and 128. 

In the table above, these values therefore indicate that:
- This read is paired 
- It mapped in a proper pair
- It's mate is on the reverse strand
- This is the second read in the pair.

You don't need to remember these codes or to work out which values add up to a specific FLAG value. 
Use the  [SAM Flag decoder](https://broadinstitute.github.io/picard/explain-flags.html). 

#### Field 6 - CIGAR Strings 
CIGAR strings describe how the read aligned to the reference. 

- `M` - alignment match (can be either a sequence match or mismatch)
-  = - Sequence match 
- `X`- Sequence mismatch
- `I` - Insertion to the reference
- `D` - Deletion from the reference
- `S` - Soft clipping

For example, the CIGAR `8M2I4M1D3M` indicates:
- 8 alignment matches
- 2 insertions
- 4 alignment matches
- 1 deletion from the reference
- 3 alignment matches
### Example

```
@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
```

The figure below shows an example of some small read alignments in **A** and what the SAM file describing these alignments would look like in **B.** 

![](https://zyxue.github.io/assets/sam_format_example.jpg) 

# VCF

- [VCF Official Documentation](https://samtools.github.io/hts-specs/VCFv4.5.pdf) 
- [Simple explanation of VCF from NYU](https://learn.gencore.bio.nyu.edu/ngs-file-formats/vcf-format/)

VCF stands for Variant Calling Format. VCF files describe variants including single nucleotide polymorphisms (SNPs), insertions and deletions (indels), and structural variants (variants bigger than 50bps). They may additionally contain genotype information for multiple samples for each variant. 

VCF files have the `.vcf` extension. 

## Format

VCF files begin with meta-information lines (lines starting with `##`) followed by a header containing the column names (starts with `#`), and finally the main body of the file. 

### Meta-information
The meta-information section contains definitions and descriptions for the data in the main body of the file. It's essential to understanding the VCF. 

Meta-information may be unstructured or structured. 
Unstructured meta-information has the format `##key=value`. In the example below, `##fileformat=VCFv4.2` and `##fileDate=20090805` are of this type, as well as a number of other lines. 

In the structured meta-information lines, the format is similar except that the value consists of a comma-separated list of key=value pairs and the entire list is enclosed within '<' and '>'. 
See the Example section for an example of how to interpret structured meta-information. 

In the main body of the VCF, there are 8 mandatory columns. 

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

The VCF below contains five different types of variants with genotype information for three samples (NA00001, NA00002 and NA00003). 

```
##fileformat=VCFv4.5
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of      Samples With    Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total  Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele   Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral       Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP     membership,     build   129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2   membership">
##FILTER=<ID=q10,Description="Quality   below   10">
##FILTER=<ID=s50,Description="Less      than    50%     of      samples have    data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype     Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype    Quality">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
20      14370   rs6054257       G       A       29      PASS    NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ     0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
20      17330   .       T       A       3       q10     NS=3;DP=11;AF=0.017     GT:GQ:DP:HQ     0|0:49:3:58,50  0|1:3:5:65,3    0/0:41:3
20      1110696 rs6040355       A       G,T     67      PASS    NS=2;DP=10;AF=0.333,0.667;AA=T;DB       GT:GQ:DP:HQ     1|2:21:6:23,27  2|1:2:0:18,2    2/2:35:4
20      1230237 .       T       .       47      PASS    NS=3;DP=13;AA=T GT:GQ:DP:HQ     0|0:54:7:56,60  0|0:48:4:51,51  0/0:61:2
20      1234567 microsat1       GTC     G,GTCT  50      PASS    NS=3;DP=9;AA=G  GT:GQ:DP        0/1:35:4        0/2:17:2        1/1:40:3
```

#### Interpreting IDs
The `INFO`and `FORMAT` columns use an ID system to store information. The meaning of each ID is stored in the meta-information at the beginning of the VCF.  IDs referenced in the `INFO` column have meta-information lines beginning with `##INFO=` while IDs referenced in the `FORMAT` column begin with `##FORMAT=`. These are both examples of structured meta-information. 

#### Interpreting structured meta-information

Let's interpret the meta-information entry below: 

```
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
```

This indicates:
- `##INFO` <- this information is relevant to the the INFO column (column 8)
- `Number=1` <- the ??????
- `ID=NS` <- The identifier for this item in this field is ``NS
- `Type=Integer` <- the type of data we are expecting in this field (a whole number in this case)
- `Description="Number of Samples with Data` <- a description of the field 

If we look at column 8 in the example, we see `NS=3` in the first entry indicating that three samples have data relevant to this variant. The third entry has `NS=2`, meaning only two samples have information relating to this variant. 

#### Interpreting sample genotype information
This example contains genotype information for three samples - NA00001, NA00002, and NA00003. We can see the sample names in the header line (the line starting `# CHROM POS ID ....`) in columns 10, 11, and 12.

The `FORMAT` column contains information on how the sample columns are formatted. In this case, the `FORMAT` column contains: 
```
GT:GQ:DP:HQ
```
To find out what these IDs mean, we look through the header for lines beginning with `##FORMAT=` and then match up the IDs. We find that:

- `GT` = genotype 
- `GQ` = genotype quality 
- `DP` = Read depth
- `HQ` = Haplotype quality

Therefore, the sample columns will list the genotype of the sample, the genotype quality, the read depth, and haplotype quality separated by colons. 
The first sample (NA00001) has the following genotype information:   `0|0:48:1:51,51`
This means:
- `0|0` - homozygous for the reference allele
- `48` - Genotype quality 48
- `1`  - read depth 1
- `51,51` - haplotype quality was 51 for both haplotypes

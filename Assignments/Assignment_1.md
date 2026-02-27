# Practical Assignment 1a - Fundamentals of Bioinformatics

This assignment covers all material in Module 1 - Fundamentals of Bioinformatics. 
It is split into two parts with 30 marks for each part. 

### Submission Details

- All bash/shell scripts should be uploaded as separate `.sh` files. 
- All other answers should be provided in a single pdf. An easy way to do this is to complete your assignment in Word (or another word processor) and export/save your final document as a pdf.  Alternatively, you are welcome to use Rmarkdown, LaTeX, or something similar that can produce pdf output. 
- Please ensure that this pdf contains your name, student number, and the assignment name (Practical Assignment 1a - Fundamentals of Bioinformatics).

## Part 1 - Bash 

You will need to produce two valid bash scripts for Part 1 of this assignment.

Q1. Write a bash script to: 
- Download the gff3 file for your assigned species (see bottom of page) to your current directory from Ensembl [1 mark]
- Count how many of each feature type there is, sorted in numerical order [4 marks]
- Count how many gene features have a gene_biotype attribute of protein_coding, and how many have a gene_biotype attribute of lncRNA [4 marks]
- Export the results to a file with a name of the form `my_species_gff_features.txt` where you use your assigned species name instead of my_species [1 mark]. NB: If the species name is missing or incorrect, no marks will be awarded for this section. The script must also include code to generate one or more comment lines in the output file/table before the table with the genome-build used, (hint: `grep` your gff to find the genome build info as the header is very large in most cases)
- Include informative comments in your script to explain each step [1 mark]
- The script must also write the code used to generate the summary (counts) data and count gene features to the output file as part of the file header. [3 marks]

Q2. For the file we used in the practicals (Drosophila_melanogaster.BDGP6.ncrna.fa), add to the final practical script provided so that:
- Two separate output files are now created, titled `dmel_ncrna_summary.txt` and `dmel_ncrna_chromosome_summary.txt` [1 mark]
- For the first output file, ensure that:
  - This output contains a meaningful header [1 mark] 
  - This output contains column names [1 marks]
  - This output includes: a) gene id; b) chromosome; c) start; d) stop; e) strand and f) gene_biotype [3 marks]
- For the second output file, ensure:
  - This output contains a meaningful header [1 mark] 
  - Count the number of gene transcripts assigned to each chromosome [4 marks]
  - Print the ID names for ncRNAs found on the Y chromosome [4]
- Include appropriate comments which make the script easier to understand [1 mark]


Note: If identical comments are identified in any submissions, a mark of **zero** will be given for this question for all suspicious submissions. 
We also strongly, **strongly** suggest completing these scripts without the use of generative AI tools. You will be able to write your assignment scripts by repurposing the code and bioinformatics tools we went over in the first three practicals. 

## Part 2 - Next Generation Sequencing

#### Data
Athaliana reference genome:
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

SRR5882797_10M_1.fastq.gz https://adelaideuniversity.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz

SRR5882797_10M_2.fastq.gz https://adelaideuniversity.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz
### 2.1. Variant Calling Script [18 marks]

 Write a bash script to perform the variant calling analysis tasks listed below. **This script should use tools and techniques that were covered in the module.** 
 
 **Note:** [1 mark] will be awarded for placing all files in the directories specified in the instructions below and [2 marks] will be awarded for including clear and concise comments throughout. 
 
  Your script should do the following:
1. Download the reference genome (i.e. fasta file) of the model plant _Arabidopsis thaliana_ to the subdirectory `ref` from the Ensembl ftp directory using the link provided above and either `curl` or `wget`. [1 mark]
2. Download the sequencing data (Illumina paired-end reads) at the links provided to the directory `0_raw` using `curl` or `wget` and rename these files to the names provided alongside the links.  [3 marks]
3. Index the reference genome for use with `bwa`. You will need to unzip it first. [1 mark]
4. Run `fastqc` on raw reads. [1 mark]
5. Trim reads for poor quality bases using `fastp`. Write trimmed reads to the directory `1_trim`. [2 marks]
6. Align paired-end reads to the genome using `bwa mem`, resulting in a single BAM file in the directory `2_align`. Don't include read groups (exclude the `-R` option and the string following it). Ensure there are no intermediary SAM (`.sam`) files saved. [3 marks]
7. Sort and index your BAM file. Remove the unsorted BAM file after you've done this. [3 marks]
8. Call variants and write the vcf to `3_variants`. You do not need to remove duplicates from your sorted BAM file and should instead use your sorted BAM file from step 7 as input for variant calling. [1 mark]
9. Filter variants to keep only variants with QUAL>=30 and print the number of variants that meet this criteria to standard output (`stdout`). [2 marks]


### 2.B Short Answer Questions [12 marks]
1. Describe the 4 main components (the four lines) of a FASTQ read/record. [2 marks]
2. Illumina short reads suffer from a deterioration in quality towards the 3' end. Describe the process which causes this. [1 mark]
3. What does it mean if reads align as a proper pair? [1 mark]
4. Illumina short reads may contain portions of adapter sequences at their 3' end. Describe how and why some reads may contain parts of an adapter while others may not. [1 marks]
5. What is the difference between a SAM and a BAM file? [1 mark] 
6. Index files are regularly encountered in bioinformatics. For example, `.bai` (and `.csi`) is an index file for BAM files and `.fai` are index files of FASTA files. Describe, in general terms, what index files facilitate. [1 mark] 
7. What is the difference between a local and global sequence alignment and under what general condition would you use them? [2 marks]
8. In protein residue sequence alignment we use a substitution matrix to score mismatches. Why can't all amino acid substitutions be considered equal? [2 marks]  
9. What does the bitwise SAM FLAG 83 indicate? [1 mark] 


## Species For Question 1

*If your student number is not listed, please contact Anna/Alex to ensure you are added to the list*

You can download your assigned species here: 'http://ftp.ensembl.org/pub/release-100/gff3/' of course you will have to add the relevant additional information to specify your species and the '.100.gff3.gz' file. 

| Student Number | Species                         | Taxonomy ID | Common Name                    |
|:---------|:--------------------------------|------------:|:-------------------------------|
|1608942|Geospiza fortis|48883|Medium Ground-Finch|
|1662634|Pan paniscus|9597|Pygmy Chimpanzee|
|1752480|Hippocampus comes|109280|Tiger Tail Seahorse|
|1774160|Equus asinus asinus|9793|Wild Ass|
|1774364|Xiphophorus maculatus|8083|Southern Platyfish|
|1822446|Stachyris ruficeps|181631|Rufous-Capped Babbler|
|1822872|Equus caballus|9796|Horse|
|1823474|Clupea harengus|7950|Atlantic Herring|
|1825565|Ovis aries|9940|Sheep|
|1828996|Gambusia affinis|33528|Western Mosquitofish|
|1851925|Gasterosteus aculeatus|69293|Three-Spined Stickleback|
|1852976|Zonotrichia albicollis|44394|White-Throated Sparrow|
|1860873|Sarcophilus harrisii|9305|Tasmanian Devil|
|1868666|Macaca fascicularis|9541|Crab-Eating Macaque|
|1871005|Gorilla gorilla|9593|Western Gorilla|
|1872208|Poecilia formosa|48698|Amazon Molly|
|1876844|Mustela putorius furo|9668|Domestic Ferret|
|1888122|Astyanax mexicanus pachon|7994|Blind Cave Fish|
|1889991|Bison bison bison|9901|Bison|
|1895247|Ictidomys tridecemlineatus|43179|Thirteen-Lined Ground Squirrel|
|1897275|Chinchilla lanigera|34839|Long-Tailed Chinchilla|
|1897522|Otolemur garnettii|30611|Small-Eared Galago|
|1899225|Cebus capucinus|9516|White-Faced Sapajou|
|1901867|Terrapene carolina triunguis|158814|Three-Toed Box Turtle|
|1902444|Mus pahari|10093|Shrew Mouse|
|1902999|Neovison vison|452646|American Mink|
|1903474|Serinus canaria|9135|Common Canary|
|1903493|Canis lupus familiaris|9612|Dog|
|1907116|Seriola dumerili|41447|Greater Amberjack|
|1914323|mus musculus c57bl6nj|10090|Mouse|
|1915786|Electrophorus electricus|8005|Electric Eel|
|1916630|Monopterus albus|43700|Swamp Eel|
|1918737|Pogona vitticeps|103695|Central Bearded Dragon|
|1921723|Anas platyrhynchos|8839|Mallard|
|1921893|Sus scrofa usmarc|9823|Pig|
|1923348|Mola mola|94237|Ocean Sunfish|
|1923580|Mus musculus lpj|10090|Mouse|
|1923733|Rattus norvegicus|10116|Norway Rat|
|1925587|Poecilia reticulata|8081|Guppy|
|1926036|Oryzias javanicus|123683|Javanese Ricefish|
|1926397|Salmo salar|8030|Atlantic Salmon|
|1934359|Anabas testudineus|64144|Climbing Perch|
|1940801|Vulpes vulpes|9627|Red Fox|
|1943794|Notechis scutatus|8663|Mainland Tiger Snake|
|1953052|Rhinolophus ferrumequinum|59479|Greater Horseshoe Bat|
|1961590|Chrysemys picta bellii|8479|Western Painted Turtle|
|1963206|Neolamprologus brichardi|32507|Cichlid|
|1963286|Castor canadensis|51338|American Beaver|
|1964213|Lepidothrix coronata|321398|Blue-Crowned Manakin|
|1976953|panthera tigris altaica|9694|Amur Tiger|
|1980273|Pelodiscus sinensis|13735|Chinese Soft-Shelled Turtle|
|1980973|Sphenodon punctatus|8508|Tuatara|
|1981829|Strigops habroptila|2489341|Kakapo|
|1983032|Vombatus ursinus|29139|Common Wombat|
|1983057|Aotus nancymaae|37293|Mas Night Monkey|
|1983060|Parambassis ranga|210632|Indian Glassy Fish|
|1983069|sus scrofa berkshire|9823|Pig|
|1983083|sus scrofa bamei|9823|Pig|
|1983089|heterocephalus glaber male|10181|African Mole-Rat|
|1983093|Callithrix jacchus|9483|White-Tufted-Ear Marmoset|
|1983110|Panthera pardus|9691|Leopard|
|1983157|sus scrofa pietrain|9823|Pig|
|1983613|Ictalurus punctatus|7998|Channel Catfish|
|1984281|Pygocentrus nattereri|42514|Red-Bellied Piranha|
|1985604|Myripristis murdjan|586833|Pinecone Soldierfish|
|1986119|Erinaceus europaeus|9365|Western European Hedgehog|
|1987833|Mus musculus akrj|10090|Mouse|
|1987856|Spermophilus dauricus|99837|Daurian Ground Squirrel|
|1988768|Monodelphis domestica|13616|Gray Short-Tailed Opossum|
|1989620|canis lupus dingo|9612|Dingo|
|1989666|Fukomys damarensis|885580|Damara Mole-Rat|
|1991177|Fundulus heteroclitus|8078|Mummichog|
|1991610|Larimichthys crocea|215358|Large Yellow Croaker|
|1991642|Poecilia mexicana|48701|Shortfin Molly|
|1992287|Gadus morhua|8049|Atlantic Cod|
|1992421|Physeter catodon|9755|Sperm Whale|
|3188346|Xiphophorus couchianus|32473|Platyfish|
|3188465|Tetraodon nigroviridis|99883|Pufferfish  |
|3191299|Tupaia belangeri|37347|Northern Treeshrew|
|3192119|Varanus komodoensis|61221|Komodo Dragon|
|3192385|Takifugu rubripes|31033|Pufferfish (fugu)|
|3198641|Taeniopygia guttata|59729|Zebra Finch|
|3198989|Pongo abelii|9601|Orangutan|
|3204028|Petromyzon marinus|7757|Lamprey|
|1983788|Junco hyemalis|40217|Dark-eyed Junco|

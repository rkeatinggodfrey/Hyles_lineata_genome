# Genome annotation from PacBio Sequel IIe assembled genome

This script describes the steps taken to annotate the Hyles lineata genome.
General assembly steps are as follows:
+ (1) Repeat modeling using Repeat Modeler http://www.repeatmasker.org/RepeatModeler/
+ (2) Repeat masking (may be redundant if you're going to run Maker annotation)
+ (3) Prep for BRAKER2 annotation: NCBI SRA toolkit and NCBI Entrez Programming Utilities t downloading trascriptomes and protein sequences from closely related species. 
+ (4) BRKAER2 annotation 

Each step includes resources used to build analysis, notes or tips on gathering data for analyses, and submission scripts

## (1) Repeat Modeler
Resources:
+ https://www.clementgoubert.com/post/a-simple-pipeline-for-te-annotation-in-an-assembled-genome
+ https://biohpc.cornell.edu/doc/annotation_2019_exercises1_v2.html
+ http://avrilomics.blogspot.com/2015/02/finding-repeats-using-repeatmodeler.html

```bash
#!/bin/bash
#SBATCH --job-name=H_lineata_repeatmod
#SBATCH -o %A_%a.220616_Hl_repeatmod.out
#SBATCH --mail-user=rkeating.godfrey@hpg.rc.ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 20
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5-00:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

date;hostname;pwd

module load repeatmodeler/2.0
module load seqkit/2.0.0

BuildDatabase -name H_lineata /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_hifiasm_220728_purge.fasta

RepeatModeler -database H_lineata -pa 32 -LTRStruct >& run.adp.out

# you can break the known and unknown repeats into separate fasta files

cat H_lineata-families.fa | seqkit fx2tab | awk '{ print "Hlineata_1.0_"$0 }' | seqkit tab2fx > H_lineata-families.prefix.fa
cat H_lineata-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > H_lineata-families.prefix.fa.known
cat H_lineata-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > H_lineata-families.prefix.fa.unknown

```


## (2) Repeat Masking

### (a) Option 1: Mask all repeats identified by Repeat Modeler

```bash
#!/bin/bash
#SBATCH --job-name=rmask_old_H_lineata
#SBATCH -o %A_%a.220728_Hl_repeatmask_old.out
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 8
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 90:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

date;hostname;pwd

module load repeatmasker/4.1.1

RepeatMasker -pa 8 -a -s -xsmall -gff -no_is -lib H_lineata-families.fa /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_hifiasm_220728_purge.fasta &> RMask.run.adp.out
```

Check and see how many repeats were masked in this single-step masking process to compare with heirarchical one that follows
```
cd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_RepeatModeler/masked_purged_v1
cat H_lineata_hifiasm_220728_purge.fasta.masked | grep -v ">" | tr -dc a-z | wc -c
```
Output: 169898824



### (b) Option 2: Mask according to suggestions from Dr. Darren Card (adapted from Yi-Ming Weng)
+ step 1: mask simple repeats
+ step 2: mask previously identified Lepidoptera-specific repeats
+ step 3: mask known repeats identified from Repeat Modeler

(1) Mask simple repeats

(2) Mask Lep repeats

(3) Mask known repeats

Check and see how many repeats were masked in this heirachical masking process as compared with masking using RepeatModeler output alone:

After step 1: 

cd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_RepeatModeler/Hl_repeatmasker_s1
```
cat H_lineata_hifiasm_220728_purge.fasta.masked | grep -v ">" | tr -dc a-z | wc -c
```
output:

After step 2:
cd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_RepeatModeler/Hl_repeatmasker_s2
```
cat H_lineata_hifiasm_220728_purge.fasta.masked.masked | grep -v ">" | tr -dc a-z | wc -c
```
output: 7222467

After step 3 (total): 
```
cd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_RepeatModeler/Hl_repeatmasker_s3
cat H_lineata_hifiasm_220728_purge.fasta.masked.masked.masked | grep -v ">" | tr -dc a-z | wc -c
```
output: 
141040752


Once satisfied with this, I renamed the final version
```mv H_lineata_hifiasm_220728_purge.fasta.masked.masked.masked H_lineata_assembly_final_3masked.fasta```

And then moved it to the parent folder
```mv H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/```

# Genome Annotation: BRAKER 2

BRAKER2 provides a method for annotating using protein or RNA-seq data. They describe using proteins from any evolutionary distance and with proteins of short evolutionary distance, but note that the latter is a decpreciated pipeline. Therefore, I used the first pipeline for my protein evidence pipeline.

Resources:
+ Running BRAKER with proteins of any evolutionary distance: 
+ Running BRAKER with RNA-seq data: https://github.com/Gaius-Augustus/BRAKER#braker-with-rna-seq-data


## (1) Running BRAKER with Protein data

## (a) Retrieve protein sequences

### First retreive protein sequences from the BUSCO arthropod database

```cd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2```

```wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
tar xvf odb10_arthropoda_fasta.tar.gz
arthropoda/Rawdata/* > arthropod.proteins.fasta```

if the wget command throws a certificate error, use:
```wget --no-check-certificate https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz```

### Retrieve protein sequences from a well-annotated, closely related species (optional) 

I additionaly downloaded M. sexta protein sequences into my ncbi downloads folder and move them into the folder where I am running BRAKER2

Download from ncbi:
```bash 
#!/bin/bash
#SBATCH --job-name=ncbi_download
#SBATCH -o %A_%a.220701_NCBI_Download.out
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 1
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 02:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

module load edirect/12.2

## Retrieve Manduca sexta protein and mRNA

esearch -db protein -query "Manduca sexta [ORGN]" | efetch -format fasta > M_sexta_protein.fasta
```

move into BRAKER2 folder:

```mv M_sexta_protein.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/ ```

Now combine these with the other arthropod proteins

```cat arthropod.proteins.fasta M_sexta_protein.fasta > all.proteins.fasta ```

## (b) Run [ProtHint](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) to create protein gff file.
BRAKER can use create a protein gff file as the first step, but it may be better to create this file on your own and then feed it into BRAKER.

```sbatch -J Hl_ProtHint prothint.sh /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/all_proteins.fasta```

```bash
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%xc_%j.log
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

dates;hostname;pwd

genome=${1}
protein=${2}

module load prothint/2.6.0

prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} ${protein}
```


```bash
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%xc_%j.log
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

dates;hostname;pwd

module load prothint/2.6.0
module load genemark_es/4.69

prothint.py /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/arthropod.proteins.fasta
```

This will create an output that is ready to use in BRAKER and AUGUSTUS:
+ ```prothint_augustus.gff```

## (c) Run BRAKER2 with protein evidence

First I ran Braker2 with protein evidence from arthropoda

When trying to run this I noticed that in the file path /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config/species there is a "Sp_1", not Hyles_lineata so I changed the name of this directory 

```sbatch -J Hl_braker2_protein Hl_braker2_protein.sh /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/prothint_augustus.gff Hyles_lineata

sbatch -J Mr_braker2_protein braker2_protein.sh /blue/kawahara/rkeating.godfrey/Manduca_rustica_genome/M_rustica_final_assembly_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/prothint_augustus.gff Manduca_rustica

```bash
#!/bin/bash
#SBATCH --job-name=%j_Hl_braker2_prot
#SBATCH --output=%j_Hl_braker2_prot.log
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
protein_gff=${2}
species=${3}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config \
--genome=${genome} --species ${species} --hints=${protein_gff} --softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```

I moved the files associated with this protein-based annoation to a folder called braker_protein_arth

## (2) Running BRAKER with RNA-seq data

### (a) Download transcriptome data as paired end reads

The ```--split-files```


```bash

#!/bin/bash
#SBATCH --job-name=H_euphorbiae_RNA
#SBATCH -o %A_%a.220616_Heuphorb_transcriptome.out
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 2
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 01:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

module load gcc/5.2.0
module load ncbi-vdb/2.8.2
module load sra/3.0.0

prefetch SRR1695429 --max-size 2500000000

# convert to FASTQ: fastq-dump will convert SRR1695429.sra to SRR1695429.fastq 
# split-files splits paired end read files 

fastq-dump --split-files  SRR1695429
```

I did this in an ncbi folder I use for downloading, so I moved the paired-end fastq files to my braker folder:

```mv *.fastq ../Hl_braker2/```

### (b) Map RNA sequencing reads to masked genome using guidance from [Kim and Kim 2022](https://www.sciencedirect.com/science/article/pii/S2666166722003860?via%3Dihub#sec2)

### if your transcriptome is a fasta file, you can specify by using -2 -f before the file path to the transcriptome file

```
#!/bin/bash
#SBATCH --job-name=%x_%j_Hl_mapreads
#SBATCH --output=%A_%a.Hl_mapreads.out
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
dates;hostname;pwd

module load hisat2
module load samtools

# create a masked genome index

hisat2-build /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta Hl_Hifi

# map RNA sequencing reads to the masked genome 

hisat2 -x Hl_Hifi -p 10 -1 /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/SRR1695429_1.fastq  \
-2 /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/SRR1695429_2.fastq | samtools sort -@ 10 -O BAM -o Hl_He_aln.bam

```

### (c) Run braker with RNA evidence

I created a new folder for this annotation evidence called braker_RNA_He and put the files needed in it (except the assembly, which lives in the parent genome folder)

```sbatch -J Hl_braker2_RNA Hl_braker2_RNA.sh /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta Hyles_lineata_RNA```

```bash
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
species=${2}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/config \
--genome=${genome} --species ${species} \
--bam=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker/braker_RNA_He/Hl_He_sort.bam \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```



### BRAKER output

The output file  braker.gtf  

The field in gft format are: 

seqname source feature start end score strand frame transcript ID and gene ID

Resources: 
+ https://github.com/Gaius-Augustus/BRAKER#output-of-braker




## (3) Evaluate gene models produced by braker2 
### (a) from arthropod protein database using BUSCO endopterygota ortholog database (odb10_lepidoptera)

```sbatch Hl_prot_model_busco.sh```


```bash
#!/bin/bash

#SBATCH --job-name=Hl_lep_prot_genemodel_busco
#SBATCH -o Hl_lep_prot_genemodel_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_arth/braker/augustus.hints.aa \
 -o ./Hl_prot_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m protein -c 12
 ```
 
 Results from dataset endopterygota_odb10
 
 endopterygota = C:94.6%[S:86.6%,D:8.0%],F:2.0%,M:3.4%,n:2124     

 lepidoptera =  C:95.6%[S:86.5%,D:9.1%],F:1.1%,M:3.3%,n:5286  

Check how many genes

```grep ">" augustus.hints.aa | wc -l```

21323


 ### (b) from Hyles euphorbiae transcriptome using BUSCO lepidoptera or endopterygota ortholog database (odb10_lepidoptera)

 ```bash
 #!/bin/bash
#SBATCH --job-name=Hl_lep_HeRNA_genemodel_busco
#SBATCH -o Hl_lep_HeRNA_genemodel_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/augustus.hints.aa \
 -o ./Hl_RNA_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m protein -c 12
 ```

Result

endopterygota = C:92.9%[S:89.2%,D:3.7%],F:3.0%,M:4.1%,n:2124 

lepidoptera = C:93.3%[S:89.4%,D:3.9%],F:1.8%,M:4.9%,n:5286 

Check how manay of genes

```grep ">" augustus.hints.aa | wc -l```

19425

# Genome Annotation: TSEBRA

### (a) Combine gene models from protein and transcriptome evidence

Resources:
+ (TSEBRA Github)[https://onlinelibrary.wiley.com/doi/epdf/10.1002/jmor.21510]

Clone the TSEBRA directory
```git clone https://github.com/Gaius-Augustus/TSEBRA```

Gather necessary files
+ augustus.hints.gtf files from Braker2 protein and RNA
+ a configuration (.cfg) file as descirbed (here)[https://github.com/Gaius-Augustus/TSEBRA#configuration-file]. You can use the default file located in /TSEBRA/config/default.ctg as a start. The default is quite strict. 


```bash
#!/bin/bash
#SBATCH --job-name=Hl_TSEBRA
#SBATCH -o Hl_TSEBRA.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 24:00:00
#SBATCH -c 24

module load python3

/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/TSEBRA/bin/tsebra.py \
--keep_gtf /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_arth/braker/augustus.hints.gtf,/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/augustus.hints.gtf \
-c /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/TSEBRA/config/default.cfg \
-e /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_arth/braker/hintsfile.gff,/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/hintsfile.gff \
-o Hl_protein_rnaseq_combine.gtf

/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/scripts/gtf2aa.pl \
/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta \
/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/Hl_protein_rnaseq_combine.gtf \
Hl_braker_final_aa.fa
```

Check how many genes:
```grep ">" Hl_braker_final_aa.fa | wc -l```

33756 this looks like too many


### (b) Run BUSCO on final gene model set

``bash
 #!/bin/bash
#SBATCH --job-name=Hl_lep_all_genemodel_busco
#SBATCH -o Hl_lep_allc_genemodel_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/Hl_braker_final_aa.fa \
 -o ./Hl_all_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12
 ```
 
 lepidoptera = C:97.8%[S:55.2%,D:42.6%],F:0.5%,M:1.7%,n:5286 

 This duplication percent seems to indicate there are a many isoforms of certain genes in this gene set. 




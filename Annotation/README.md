# Genome annotation from PacBio Sequel IIe assembled genome

This script describes the steps taken to annotate the Hyles lineata genome.
General assembly steps are as follows:
+ (1) Repeat modeling using Repeat Modeler http://www.repeatmasker.org/RepeatModeler/
+ (2) Repeat masking
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
+ Output files: https://github.com/Gaius-Augustus/BRAKER#output-of-braker

Braker will create a subfolder called "braker" within which you fill find:
+ braker.gtf with the field names seqname source feature start end score strand frame transcript ID and gene ID
+ augustus.hints.aa file that contains transcripts. Run BUSCO on these using the Lepidoptera database.




## (1) Running BRAKER with Protein data

## (a) Retrieve protein sequences

### Retrieve protein sequences from a well-annotated, closely related species

I downloaded M. sexta protein sequences into my ncbi downloads folder and move them into the folder where I am running BRAKER2

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

check how many proteins are in this file: ```grep ">" M_sexta_protein.fasta | wc -l```

53159

move into BRAKER2 folder:

```mv M_sexta_protein.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/ ```


## (b) Run [ProtHint](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) to create protein gff file.
BRAKER can use create a protein gff file as the first step, but it may be better to create this file on your own and then feed it into BRAKER.

```sbatch -J Hl_ProtHint prothint.sh /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/M_sexta_protein.fasta```

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
module load genemark_es/4.69

prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} ${protein}
```

I ended up submitting the prothint submission script below instead. I got an error for not having genemark installed

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

prothint.py /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/M_sexta_protein.fasta
```

This will create an output that is ready to use in BRAKER and AUGUSTUS:
+ ```prothint_augustus.gff```

## (c) Run BRAKER2 with protein evidence

First I ran Braker2 with protein evidence from Manduca sexta

When trying to run this I noticed that in the file path /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config/species there is a "Sp_1", not Hyles_lineata so I changed the name of this directory 

```sbatch -J Hl_prot_braker2 Hl_braker2_protein.sh /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_Ms/prothint_augustus.gff Hyles_lineata```

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



## (2) Running BRAKER with RNA-seq data

### (a) Map RNA-Seq data from caterpillar and adult to genome

(RNA reads from Jay Goldberg, University of Arizona)

First you need to trim using trimmomatic. If reads are from ncbi they may already be trimmed.

sbatch -J Hl_trimmomatic trimmomatic.sh

```bash
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load trimmomatic/0.39

for sample in $(ls *fq.gz | cut -d "_" -f 1,2 | sort | uniq) # the -f in this line indicates "retain fields 1 and 2 deliminted by an underscore
do
    fq1=$(ls ${sample}_1*)
    fq2=$(ls ${sample}_2*)
    trimmomatic PE -threads 16 \
    ${fq1} ${fq2} \
    ${sample}_1_clean.fq.gz ${sample}_1_unpaired.fq.gz \
    ${sample}_2_clean.fq.gz ${sample}_2_unpaired.fq.gz LEADING:3 TRAILING:3 MINLEN:36
done
```


Following trimming, I put all of the trimmed filed in a subfolder called "trimmed" (not necessary if you do not want to). You will map these reads to your genome using Hisat2 and then convert the output sam files to bam files using samtools.

Resources:
+ [RNA-seq analysis in R](https://bioinformatics-core-shared-training.github.io/RNAseq_September_2019/html/C_Alignment_with_HISAT2_practical.html)


```sbatch -J Hl.hisat hisat2.sh```

```bash
#!/bin/bash
#SBATCH --job-name=%x_%j_Hl_mapreads
#SBATCH --output=%A_%a.log
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

module load hisat2/2.2.1-3n
module load samtools 

hisat2-build /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/trimmed/Hlhisat

for sample in $(ls /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/trimmed/*fq.gz | cut -d "_" -f 1,2,3,4,5,6,7 | sort | uniq)
do
    fq1=$(ls ${sample}_1_clean*)
    fq2=$(ls ${sample}_2_clean*)
    name=$(echo ${sample} | cut -d "_" -f 1,2)
    hisat2 -p 32 \
    -x /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/Hlhisat \
    -1 ${fq1} -2 ${fq2} --phred33 --rna-strandness FR | samtools sort -@ 10 -O BAM -o Hl_Hl_aln.bam
done

```
The .log file output will contain mapping percentages. For the RNA-seq files the mapping rates were as follows:
+ Adult Male: 84.50%
+ Adult Female: 87.80%
+ Larva: 91.36%


### (b) Run BRAKER with RNA-seq bam alignment file

```sbatch -J Hl_braker2_RNAseq braker2.RNAseq.sh /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta Hyles_lineata_RNAseq```

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
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config \
--genome=${genome} --species ${species} \
--bam=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/Hl_Hl.aln.bam \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```


## (3) Running BRAKER with assembled transcriptome from closely related species (in this case Hyles euphorbiae)

### (a) Download transcriptome data as paired end reads

The ```--split-files``` in this script provides paired end reads

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
-2 /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/SRR1695429_2.fastq | samtools sort -@ 10 -O BAM -o Hl_He_sort.bam

```

### (c) Run braker with transcriptome evidence

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
#SBATCH --cpus-per-task=24
dates;hostname;pwd

genome=${1}
species=${2}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config \
--genome=${genome} --species ${species} \
--bam=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/Hl_He_sort.bam \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```



## (3) Evaluate gene models produced by braker2 using BUSCO Lepidoptera ortholog database (odb10_lepidoptera)

### (a) from Manduca sexta protein database 

```sbatch Hl_Ms_prot_model_busco.sh```


```bash
#!/bin/bash
#SBATCH --job-name=Hl_lep_all_genemodel_busco
#SBATCH -o Hl_lep_all_genemodel_busco.log
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
busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_Ms/braker/augustus.hints.aa \
 -o ./Hl_Ms_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12
 ```
 
Result:

lepidoptera = C:95.9%[S:85.8%,D:10.1%],F:0.9%,M:3.2%,n:5286 

How many genes predicted using Manduca sexta protein evidence?
```grep ">" augustus.hints.aa | wc -l```

20967

I moved the files associated with this protein-based annoation to a folder called braker_prot_Ms



 ### (b) from Hyles lineata RNAseq evidence 

```sbatch Hl_RNAseq_model_busco.sh```

```bash
#!/bin/bash
#SBATCH --job-name=Hl_lep_HlRNA_genemodel_busco
#SBATCH -o Hl_lep_HlRNA_genemodel_busco.log
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
busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/braker/augustus.hints.aa \
 -o ./Hl_RNA_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12 
 ```

Result:

lepidoptera =  C:95.7%[S:84.3%,D:11.4%],F:1.4%,M:2.9%,n:5286 

How many genes predicted using Hyles lineata RNA-seq evidence?
```grep ">" augustus.hints.aa | wc -l```

20268

I moved the files associated with this transcriptome-based annoation to a folder called braker_RNA_Hl


 
 ### (c) from Hyles euphorbiae transcriptome evidence

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

Result:

lepidoptera = C:93.3%[S:89.4%,D:3.9%],F:1.8%,M:4.9%,n:5286 

How many genes predicted using Hyles euphorbiae transcriptome as evidence?

```grep ">" augustus.hints.aa | wc -l```

19425

I moved the files associated with this transcriptome-based annoation to a folder called braker_RNA_He



# Genome Annotation: TSEBRA

TSEBRA is a transcript selector that allows you to combine the Braker outputs from different evidence into a single list of transcripts / predicted amino acid sequences.

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
--keep_gtf /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_Ms/braker/augustus.hints.gtf,/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/braker/augustus.hints.gtf,/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/augustus.hints.gtf \
-c /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/TSEBRA/config/default.cfg \
-e /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_prot_Ms/braker/hintsfile.gff,/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/braker/hintsfile.gff,/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/hintsfile.gff \
-o Hl_all_combine.gtf

/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/scripts/gtf2aa.pl \
/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_assembly_final_3masked.fasta \
/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/Hl_all_combine.gtf \
Hl_all_tsebra_aa.fa
```

Check how many genes:
```grep ">" Hl_all_tsebra_aa.fa | wc -l```

18336 (with default.ctg set)

This number of genes seems low, so I am going to try to run TSEBRA with a different config file, namely one that prioritizes the RNA-seq evidence.
The cfg file looks like this:

```# Weight for each hint source
# Values have to be >= 0
P 0.1
E 10000
C 5
M 1
# Required fraction of supported introns or supported start/stop-codons for a transcript
# Values have to be in [0,1]
intron_support 0.25
stasto_support 2
# Allowed difference for each feature 
# Values have to be in [0,1]
e_1 0.25
e_2 1
# Values have to be >0
e_3 25
e_4 10
```

Check how many genes in this outup
```grep ">" Hl_all_tsebra_prefb1_aa.fa | wc -l```

20003 (when RNA evidence is prefered)


### (b) Run BUSCO on this gene model set

```bash
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
busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_He/braker/Hl_all_tsebra_prefb1_aa.fa \
 -o ./Hl_all_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12
 ```
 
Result:

lepidoptera = C:95.9%[S:80.7%,D:15.2%],F:0.6%,M:3.5%,n:5286  

This seems like a pretty good result to me. The duplication is a little high, but I will move forward with this as my final predicted transcript set from structural annotation. The next step will be to functionally annotate these transcripts / putative genes. 

UPDATE: Braker may be predicting too many transcripts from related organism protein evidence or TSEBRA may not be properly culling duplicate transcripts and therefore I will use only RNAseq-based annotation for functional steps and color gene analysis 

# Functional Annotation: DIAMOND

Diamond is a sequences aligner that uses blastp to search for your translated transcripts (braker's aa output) in a database of sequences (chosen by you). We will blast our transcripts against two databases: (1) ncbi's non-redundant protein database and (2) the Swiss Prot arthropod database.

Resources:
+ https://github.com/bbuchfink/diamond
+ https://github.com/bbuchfink/diamond/wiki
+ tutorial: https://github.com/bbuchfink/diamond/wiki/1.-Tutorial


### (1) Non-redundant protein database
you can download the database from ncbi here
```wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"```
I used an existing copy on the orange (storage) drive of the UF hipergator cluster (provided by postdoc YiMing Weng)

You will also turn the output .tsv files into .gff files using the python script blast2gff.py available from [genomeGTFtools](https://github.com/wrf/genomeGTFtools#blast2gff). Citation [here](https://elifesciences.org/articles/31176). 

```sbatch -J Hl_nr /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/diamond/Hl_diamond.sh /orange/kawahara/yimingweng/databases/nr.dmnd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/braker_RNA_Hl/braker/augustus.hints.aa 0.00001 Hl_RNA_nr_k5_1e5```

```bash
#!/bin/bash
#SBATCH --job-name=Hl_diamond
#SBATCH -o Hl_diamond.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 48:00:00
#SBATCH -c 24
#SBATHC --account=kawahara
#SBATCH --qos=kawahara-b

module load diamond/2.0.9
module load python3

# example

database=${1} # full path to the database in dmnd format, the diamond will use it to find the function for the querying gene model
gene_model=${2} # gene model for functional annotation in fasta format
cutoff=${3}
outname=${4}
path=$(echo ${database} | rev | cut -d "/" -f 2- | rev)
database_name=$(echo ${database} | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
dmnd=$(ls ${path}/${database_name}.dmnd)

if [ -z "${dmnd}" ]
then
  echo "converting database from fasta to dmnd format"
  diamond makedb --in ${database} -d nr
else
  echo -e "the database has been converted to dmnd format, skip this step and run diamond"
  diamond blastp -k5 -e ${cutoff} -d ${database} -q ${gene_model} -o ${outname}.tsv
fi

# convert the tsv file to gff file
python3 /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/diamond/blast2gff.py -b ${outname}.tsv > ${outname}.gff
```

### (2) Arthropod uniprot database

Here I also link to an existing copy of the database on the HiperGator Orange drive

/orange/kawahara/yimingweng/databases/uniprot_arthropod.dmnd

sbatch -J Hl_uniprot /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/diamond/Hl_diamond.sh /orange/kawahara/yimingweng/databases/uniprot_arthropod.dmnd /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_braker2/Hl_all_tsebra_prefb1_aa.fa 0.00001 Hl_uniprot_k5_1e5


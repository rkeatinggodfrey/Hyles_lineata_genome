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
#SBATCH --mail-type=FAIL,END
#SBATCH -c 20
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5-00:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

date;hostname;pwd

module load repeatmodeler/2.0

BuildDatabase -name H_lineata /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_hifiasm_220728_purge.fasta

RepeatModeler -database H_lineata -pa 32 -LTRStruct >& run.adp.out

# you can break the known and unknown repeats into separate fasta files

cat H_lineata-families.fa | seqkit fx2tab | awk '{ print "Hlineata_1.0_"$0 }' | seqkit tab2fx > H_lineta-families.prefix.fa
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

RepeatMasker -pa 8 -a -s -xsmall -gff -no_is -lib H_lineata-families.fa /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_hifiasm_220728_purge.fasta$
```

### (b) Option 2: Mask according to suggestions from Dr. Darren Card (adapted from Yi-Ming Weng)
+ step 1: mask simple repeats
+ step 2: mask previously identified Lepidoptera-specific repeats
+ step 3: mask known repeats identified from Repeat Modeler

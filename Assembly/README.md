# Hyles_lineata_genome
## Assembly of the white-lined sphinx moth genome from PacBio Sequel IIe long reads

## Orientation to Pac Bio Sequel IIe data files
PacBio Sequel IIe will provide you a number of files as outputs.
A description of files can be found here: https://www.pacb.com/wp-content/uploads/Sequel_II_and_IIe_Data_Files.pdf

The primary files of interest to you are:

(1) yourjobnumber.bam # PacBio Hifi Reads (>=QV 20); additional filtering for read quality; can be done using the rq tag

(2) yourjobnumber.fastq.gz # contains sequence information with stats for each base: Same reads as .bam file but with less information about individual reads. Can be used directly in downstream applications using analysis tools for HiFi Reads.

(3) yourjobnumber.fasta.gz # contains sequences without stats/metadata

## Fast QC on reads

The fastqc on the fastq file took about 45 minutes to run, so I suggest submitting it as a SLURM job


```bash
#!/bin/bash
#SBATCH --job-name=H_lineata_fastqc
#SBATCH --output=H_lineata_raw_reads_fastqc-%j.out
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=02:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

## load fastqc module
module load fastqc

## run fastqc on the data
fastqc m64219e_220329_140935.hifi_reads.fastq.gz

```


## PacBio adapter sequence check

Source: commands provided by Yi-Ming Weng, postdoc in Kawahara lab

```zcat m64219e_220329_140935.hifi_reads.fastq.gz | grep -v "@" | grep "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"  | wc -l```

This is for the adapter, where the sequence in grep is the adapter sequence in the UniVec database of NCBI
output: ### found 395 matches

```zcat m64219e_220329_140935.hifi_reads.fastq.gz | grep -v "@" | grep "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA"  | wc -l```
This is for the Pacific Biosciences C2 Primer
output: ### found 0 match


## HiFiasm Genome Assembly 

Source: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html

This version uses the most aggressive duplicate purging using the option -l 3

```bash
#!/bin/bash
#SBATCH --job-name=H_lineata_assembly_hifiasm
#SBATCH -o %A_%a.220728_Hl_assembly.out
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 30
#SBATCH --mem-per-cpu=5gb
#SBATCH -t 30:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara-b


date;hostname;pwd

module load ufrc
module load hifiasm

hifiasm -o /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_assembly/H_lineata_hifiasm_220728.asm -l 3 -t 30 /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/m64219e_220329_140935.hifi_reads.fastq.gz
```

## Genome assembly quality assessment

Description of step:

Here we use the assemblystat.py script (source https://www.nature.com/articles/s41592-020-01056-5)

You will need to convert your primary contig file from .gfa format into .fa (fasta) format to run the assemblystat script.

Run this line of code in the terminal to produce a FASTA file from the GFA file from primary contig ctg file 
```awk '/^S/{print ">"$2;print $3}' H_lineata_hifiasm_220728.asm.bp.p_ctg.gfa > H_lineata_hifiasm_220728.asm.bp.p_ctg.fa```

If there is already a copy of the assembly.py script in another location, you can copy it to your folder using
```cp -r /file/path/to/source/location/assemblystats.py .```

Otherwise, there is a copy of the assemblystats.py script lives in this directory

Change the permissions so you can execute the file 

```chmod +x assemblystats.py ```

And load the python module to run it

```module load python```



## Genome size from kmers

Resources: 
+ https://github.com/tbenavi1/genomescope2.0
+ https://doi.org/10.1093/bioinformatics/btx304
+ https://github.com/refresh-bio/KMC

Citations:
 + Kokot et al. 2017 https://doi.org/10.1093/bioinformatics/btx304
 + Deorowicz et al. 2015 https://doi.org/10.1093/bioinformatics/btv022 Pages 1569–1576
 + Deorowicz et al. 2013 https://doi.org/10.1186/1471-2105-14-160

```bash
#!/bin/bash
#SBATCH --job-name=Hlineata_kmc
#SBATCH -o Hlineata_kmc_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH -c 3
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 00:30:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

module load kmc/3.2.1

# create directory for kmc temporary files
mkdir kmc_tmp
 
kmc -k29 /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/m64219e_220329_140935.hifi_reads.fastq.gz 29mers kmc_tmp

# Having the k-mers counted it is possible to dump KMC binary database to textual form with kmc_tools.

kmc_tools transform 29mers dump 21mers.txt

kmc_tools transform 29mers histogram 21mer_reads.histo
```

After generating the .histo file, it can be dragged into the [GenomeScope GUI](http://qb.cshl.edu/genomescope/genomescope2.0/) to produce k-mer profile



## BUSCO 

If you prefer, you can make a directory for your BUSCO submission scripts and outputs

```mkdir Gs_busco ```

Where "Gs" is an abbreivation of the genus and species

I copied the config.ini file from another folder on our cluster

```cp /orange/kawahara/amanda.markee/prelim_luna_dge/genome_assembly/spades/spades_assembly_all_illumina/K77/BUSCO/config.ini . ```

Edit SLURM submission script to ensure that right paths are set for your genome assembly and busco contig file

Make sure you load the latest version of busco or the most recent one on the cluster. You can use the command "module spider busco" to see what versions are available on the cluster

The first two lines of this script (BUSCO_CONFIG_FILE and AUGUSTUS_CONFIG_FILE) should point to the configuration file and the busco folder, respectively

The busco command should be followed the path to your assembled genome fasta file and should end with a " \" parameters for BUSCO are as follows:

 + -i input file line
 + -o output folder
 + -l lineage dataset (endopterygota_odb10, insecta_odb10, or lepidoptera_odb10)
 + -m mode set to "genome", "proteins", or "transcriptome"


```bash
#!/bin/bash
#SBATCH --job-name=Hlineata_busco
#SBATCH -o Hlineata_busco5_endo_ctg_220728%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 3:00:00
#SBATCH -c 6
#SBATCH --account=kawahara

export BUSCO_CONFIG_FILE=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/config.ini
export AUGUSTUS_CONFIG_PATH=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/

echo $BUSCO_CONFIG_FILE

module load busco/5.2.0

busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_assembly/H_lineata_hifiasm_220728.asm.bp.p_ctg.fa \
 -o BUSCO_Hlineata_endopterygota -l /data/reference/busco/v5/lineages/endopterygota_odb10        \
 -m genome -c 6
 ```

```Results: C:99.2%[S:95.9%,D:3.3%],F:0.3%,M:0.5%,n:2124```


 ## Duplicate purging

Source: Script from Yi-Ming Weng, a postdoc in the Kawahara lab following purge_haplotig pipeline (https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)

 If the BUSCO score indicate duplicated genes > 1%, you may need to purge duplicates.

 Overview of method:
 + (1) Map raw reads (subreads) to the genome assembly
 + (2) Take a look a the resulting histogram and decide the cut offs
 + (3) Purge duplicates
 + (4) Re-run BUSCO on duplicate purged assembly

### (1) Map raw reads (subreads) to the genome assembly (if you ran minimap2 for blobplot, you already have these as aligned bam files and only need to run the "purge_haplotigs" part of the code).

```bash
#!/bin/bash
#SBATCH --job-name=Hl_minimap
#SBATCH -o Hl_Kely_minimap.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

## This part of the script makes the aligned .bam files

minimap2 -t 4 -ax map-pb /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/m64219e_220329_140935.hifi_reads.fasta \
/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/m64219e_220329_140935.hifi_reads.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o Hl.aln.bam -T tmp.ali

## This part of the script generates a png of a coverage histogram and a genomecov file that you'll use in the second step 
purge_haplotigs  hist  \
-b //blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_blob/bam.mini/Hl.aln.bam  \
-g /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/m64219e_220329_140935.hifi_reads.fasta 
```

### (2) Determine cut offs

Take a look a the resulting histogram and decide the cut offs 
+ -l read depth low cutoff
+ -m low point between haploid and diploid peaks
+ -h read depth high cutoff

+ -j auto-assign contig as junk if this % or greater is low/high coverage
+ -s auto-assign contig as suspected haplotig if the % or less of the contig is diploid level of coverage

I decided 6 an 118 with a median value of 34


```bash
#!/bin/bash
#SBATCH --job-name=Hl_purge_s2
#SBATCH -o Hl_s2_purge_cut.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl.aln.sorted.bam.gencov  \
-l 6  \
-m 34  \
-h 118  \
-o coverage_stats.csv \
-j 80 \
-s 80
```

### (3) Purge duplicates

Input your original HiFiasm assembly and the coverage_stats.csv output from step 2.

```bash
#!/bin/bash
#SBATCH --job-name=Hl_s3_purge_
#SBATCH -o Hl_purge_haplotigs.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs purge  \
-g /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_assembly/H_lineata_hifiasm_220728.asm.bp.p_ctg.fa  \
-c /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/coverage_stats.csv \
-o H_lineata_hifiasm_220728_purge
```

Compare he previous assembly size with the duplicate-purged assembly 

Hyles lineata went from 471108619 to 452617722


### (4) Rerun BUSCO

Run BUSCO on duplicate-purged genome 

```bash
#!/bin/bash
#SBATCH --job-name=Hlineata_busco
#SBATCH -o Hlineata_busco5_endo_ctg_220728%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 3:00:00
#SBATCH -c 6
#SBATCH --account=kawahara

export BUSCO_CONFIG_FILE=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/config.ini
export AUGUSTUS_CONFIG_PATH=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/

echo $BUSCO_CONFIG_FILE

module load busco/5.2.0

busco -f -i /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_assembly/H_lineata_220728.purge.fasta \
 -o BUSCO_Hlineata_lepidoptera -l /data/reference/busco/v5/lineages/lepidoptera_odb10        \
 -m genome -c 6
 ```

```Results: C:98.9%[S:97.9%,D:1.0%],F:0.2%,M:0.9%,n:5286```

Looks like this worked to reduce the duplication percent while not reducing the completeness BUSCO a lot

 Check and make sure adapter sequences are not in the assembly
```cat H_lineata_hifiasm_220728_purge.fasta | grep -v "@" | grep "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" | wc -l```


## Blobtools

Run a contamination check using blobtools

Scripts and guidance from Shova Mishra (shovamishra@ufl.edu), postdoc in Peter DiGennaro lab 

Resources:
+ https://blobtools.readme.io/docs/seqfilter
+ https://github.com/DRL/blobtools
+ https://blobtools.readme.io/docs/view

### Files needed to run Blobtools
To run blob tools we need three files (b and c can be submitted as jobs at the same time)
+ (a) nodes.dmp and names.dmp files (from NCBI tax dump) 
+ (b) .bam alignment file (from minimap2 or bowtie2 script) 
+ (c) .nt blast match file (from NCBI megablast)


### (a) NCBI tax dump

To download NCBI taxdump and create nodes.dmp and names.dmp, navigate to the folder where you want the data (mine is called Hl_blob), then copy and paste this command in terminal. It will create a data directory with files inside

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/ \
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp \
./Hl_blob nodesdb --nodes data/nodes.dmp --names data/names.dmp

### (b) Mapped reads .bam file from minimap2

Resources:
+ https://github.com/lh3/minimap2#general

Note that I use the -ax parameter set to map-pb based on advice from this github for pacbio hifi reads

```bash
#!/bin/sh
#SBATCH --job-name=Hl_minimap2
#SBATCH --output=Hl_minimap2_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60gb
#SBATCH --time=08:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

pwd; hostname; date

module load minimap2

minimap2 -ax map-pb H_lineata_220728.purge.fasta m64219e_220329_140935.hifi_reads.fastq.gz > Hl.aln.sam

module load samtools
#convert SAM file to BAM file
samtools view -S -b Hl.aln.sam > Hl.aln.bam

#Use samtools sort to convert the BAM file to a coordinate sorted BAM file
samtools sort Hl.aln.bam > Hl.aln.sorted.bam

#index a genome sorted bAM file for quick alignment
samtools index Hl.aln.sorted.bam > Hl_indexed_sorted_bam
```

### (c) megablast.nt file of matched hits

While you are running minimap2, you can also run megablast. 

```bash
#!/bin/sh
#SBATCH --job-name=Hl_megablast
#SBATCH --output=Hl_megablast_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb
#SBATCH --time 8-00:00:00
#SBATCH --qos=kawahara
#SBATCH --account=kawahara

pwd; hostname; date

module load ncbi_blast
blastn -db nt -task megablast -query /blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/H_lineata_hifiasm_220728_purge.fasta -out Hlineata_megablast.nt -evalue 1e-5 -outfmt "6 qseqid staxids bitscore sgi sskingdoms sscinames" -max_target_seqs 1 -num_threads=16
```

### Now use these files to run blobtools and plot the results

```bash
#!/bin/sh
#SBATCH --job-name=Hl_blob
#SBATCH --output=Hl_blob_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rkeating.godfrey@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30gb
#SBATCH --time 04:00:00
#SBATCH --qos=kawahara
#SBATCH --account=kawahara

pwd; hostname; date

## Run blob tools
 
module load blobtools/1.0

blobtools create -i H_lineata_hifiasm_220728_purge.fasta -b Hl.aln.bam -t Hlineata_megablast.nt --nodes nodes.dmp --names names.dmp -o Hlineata
## You can then view and plot
blobtools view -i Hlineta.blobDB.json
blobtools plot -i Hlineata.blobDB.json

```


## Assembly visualization

Make some nice plots of your assembly stats
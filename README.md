# Hyles_lineata_genome
Assembly and annotation of white-lined sphinx moth genome

# Orientation to Pac Bio Sequel IIe data files
PacBio Sequel IIe will provide you a number of files as outputs.
A description of files can be found here: https://www.pacb.com/wp-content/uploads/Sequel_II_and_IIe_Data_Files.pdf

The primary files of interest to you are:

(1) yourjobnumber.bam # PacBio Hifi Reads (>=QV 20); additional filtering for read quality; can be done using the rq tag

(2) yourjobnumber.fastq.gz # contains sequence information with stats for each base: Same reads as .bam file but with less information about individual reads. Can be used directly in downstream applications using analysis tools for HiFi Reads.

(3) yourjobnumber.fasta.gz # contains sequences without stats/metadata



#############################
#### FASTQC ON SEQUENCES ####
#############################

# The fastqc on the fastq file took about 45 minutes to run, so I suggest submitting it as a SLURM job

########-----START FASTQC ASSEMBLY SCRIPT CONTENT-----######
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

########-----END FASTQC ASSEMBLY SCRIPT CONTENT-----######

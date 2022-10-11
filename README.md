# Hyles_lineata_genome
# Assembly and annotation of white-lined sphinx moth genome

## Orientation to Pac Bio Sequel IIe data files
PacBio Sequel IIe will provide you a number of files as outputs.
A description of files can be found here: https://www.pacb.com/wp-content/uploads/Sequel_II_and_IIe_Data_Files.pdf

The primary files of interest to you are:

(1) yourjobnumber.bam # PacBio Hifi Reads (>=QV 20); additional filtering for read quality; can be done using the rq tag

(2) yourjobnumber.fastq.gz # contains sequence information with stats for each base: Same reads as .bam file but with less information about individual reads. Can be used directly in downstream applications using analysis tools for HiFi Reads.

(3) yourjobnumber.fasta.gz # contains sequences without stats/metadata

# Fast QC on reads

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


# PacBio adapter sequence check

Source: commands provided by Yi-Ming Weng, postdoc in Kawahara lab

```zcat m64219e_220329_140935.hifi_reads.fastq.gz | grep -v "@" | grep "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"  | wc -l```

This is for the adapter, where the sequence in grep is the adapter sequence in the UniVec database of NCBI
output: ### found 395 matches

```zcat m64219e_220329_140935.hifi_reads.fastq.gz | grep -v "@" | grep "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA"  | wc -l```
This is for the Pacific Biosciences C2 Primer
output: ### found 0 match


# HiFiasm Genome Assembly 

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

# Genome assembly quality assessment

Description of step:

Here we use the assemblystat.py script (source https://www.nature.com/articles/s41592-020-01056-5)

You will need to convert your primary contig file from .gfa format into .fa (fasta) format to run the assemblystat script.

Run this line of code in the terminal to produce a FASTA file from the GFA file from primary contig ctg file 
```awk '/^S/{print ">"$2;print $3}' H_lineata_hifiasm_220728.asm.bp.p_ctg.gfa > H_lineata_hifiasm_220728.asm.bp.p_ctg.fa```

If there is already a copy of the assembly.py script in another location, you can copy it to your folder using
```cp -r /file/path/to/source/location/assemblystats.py .```

Change the permissions so you can execute the file 

```chmod +x assemblystats.py ```

And load the python module to run it

```module load python```


## assemblystat.py script content

```python
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
```

# BUSCO 

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


 # Duplicate purging

Source: Script from Yi-Ming Weng, a postdoc in the Kawahara lab following purge_haplotig pipeline (https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)

 If the BUSCO score indicate duplicated genes > 1%, you may need to purge duplicates.

 Overview of method:
 + (1) Map raw reads (subreads) to the genome assembly
 + (2) Take a look a the resulting histogram and decide the cut offs
 + (3) Purge duplicates
 + (4) Re-run BUSCO on duplicate purged assembly

## (1) Map raw reads (subreads) to the genome assembly (if you ran minimap2 for blobplot, you already have these as aligned bam files and only need to run the "purge_haplotigs" part of the code).

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

## (2) Determine cut offs

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

## (3) Purge duplicates

Input your original HiFiasm assembly adn the coverage_stats.csv output from step 2.

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

# Genome size from kmers


 


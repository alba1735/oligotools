# oligoTools

OligoTools is a Python package for designing and analyzing oligonucleotides for Northern blots, RT-qPCR, and other applications. It will generate a set of unique oligos in a size range and will create a bigwig file of the oligo coverage across a genome. It can also be used to analyze the coverage of a set of oligos across a genome.

## Installation

OligoTools can be installed using conda:

```bash
conda env create -f environment.yml
conda activate oligoTools
```

## Requirements

OligoTools requires the following:
* The [gtRNAdb](http://gtrnadb.ucsc.edu/) fasta and bed files for your genome of interest

Optional requirements:
* A reference genome (.fa) for your genome of interest
    * This can be downloaded from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html)
* A genome index (.fa.fai) for your genome of interest to generate a bigwig file for genome browser visualization
    * This can be generated using `samtools faidx <genome.fa>`
* A blast database to check for off-targets in your genome of interest
    * This can be generated using `makeblastdb -in <genome.fa> -dbtype nucl -out <genome>`


## Usage

### Generate oligos

This tool can be used to generate oligos from a set of files from gtRNAdb. After downloading the files from gtRNAdb, the files can be processed using the following command:

```bash
python oligoTools.py generate -f <gtRNAdb fasta file> -b <gtRNAdb bed file> -o <output directory (default: oligos)>
```

This will generate a csv and bed files containing possible oligo targets. The following flags can also be used:
* `-k` or `--kmerrange`: The size range of oligos to generate (default: 20-30)
* `-g` or `--genomeindex`: A genome index (.fa.fai) to use for generating the bigwig file for genome browser visualization
* `--log`: A log file to write the output to

### Analyze oligos

Once a series of target oligos have been selected, the coverage of these oligos across a genome can be analyzed. Place the target oligos in a `.tsv` file that looks like this:

```tsv
Ala-AGC-1-1_5p  CGCTCTAGGATATGAGCTAATCCC
iMet-CAT-1-1_5p  tRNA-iMet-CAT-1-1_target_0
Custom_target CGCTCTAGGATATGAGCTAATCCA
```

 The first column will be the oligo name for reference and ordering and the second column can contain the oligo target name from the previous step, a sequence from the previous step, or a custom sequence. The following command can be used to analyze the coverage of these oligos across a genome:

```bash
python oligoTools.py analyze -f <gtRNAdb fasta file> -l <oligo list file (all_oligos.csv)> -t <oligo target file> -o <output directory (default: oligos)>
``` 

This will generate a series of text files containing the oligos, if they map properly to the tRNA targets, where they map along the tRNA and if they are mapping to extra tRNA targets. Additionaly a `combine_IDT.csv` file is generated that can be used to quickly bulk order the oligos from IDT. 

The following flags can also be used:
* `-d` or `--idtconfig`: A config file containing the IDT order information to truncate to the end (default: "/3bio/,100nm,HPLC")
* `-b` or `--blastdb`: A blast database can also be provided to check for off-targets and this information will be added to the output files
* `-g` or `--bedgtf`: A bed file created from a gtf file to check what blast hits are mapping to
    * This can be generated using `awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$14,$7}}' <gtf file> | tr -d '";' > <gtf file>.bed` if you only want to look at genes or `awk 'OFS="\t" {print $1,$4-1,$5,$10,$14,$7}' <gtf file> | tr -d '";' > <gtf file>.bed` if you want to look at all features from the gtf file
* `--log`: A log file to write the output to
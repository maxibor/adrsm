<img src="./img/logo_adrsm.png" width="300">

# Introduction
ADRSM (Ancient DNA Read Simulator for Metagenomics) is a tool designed to simulate the paired-end sequencing of a metagenomic community. ADRSM allows you to control precisely the amount of DNA from each organism in the community, which can be used to benchmark different metagenomics methods.

# Dependencies
- [Conda](https://conda.io/miniconda.html)  

# Installation

```
conda install -c maxibor adrsm
```

# Usage

```
adrsm -d path/to/genome/directory configFile.txt
```

# Help

```
maxime@gph:~$ adrsm --help
usage: ADRSM [-h] [-d DIRECTORY] [-r READLENGTH] [-l LENSTDEV] [-fwd FWDADAPT]
             [-rev REVADAPT] [-e ERROR] [-o OUTPUT] [-s STATS]
             confFile

Ancient DNA Read Simulator for Metagenomics

positional arguments:
  confFile       path to configuration file

optional arguments:
  -h, --help     show this help message and exit
  -d DIRECTORY   path to genome directory. Default = .
  -r READLENGTH  Average read length. Default = 76
  -l LENSTDEV    Insert length standard deviation. Default = 10
  -fwd FWDADAPT  Forward adaptor. Default = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
                 NNNNNNATCTCGTATGCCGTCTTCTGCTTG
  -rev REVADAPT  Reverse adaptor. Default =
                 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
  -e ERROR       Illumina sequecing error. Default = 0.01
  -o OUTPUT      Output file basename. Default = ./metagenome.*
  -s STATS       Statistic file. Default = stats.csv

```

## Genome directory
Each genome `fasta` file must be names after the name of the organism. (ex: `Bacillus_anthracis_genome.fa`)

## Configuration file (`confFile`)
The configuration file is a `.csv` file describing, one line per genome, the mean insert size, and the expected genome coverage.
Example [short_genome_list.csv](./short_genome_list.csv):

```
genome, insert_size, coverage
Agrobacterium_tumefaciens_genome.fa, 47 , 0.1
Bacillus_anthracis_genome.fa, 48, 0.2
```

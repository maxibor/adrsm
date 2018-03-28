[![Anaconda-Server Badge](https://anaconda.org/maxibor/adrsm/badges/installer/conda.svg)](https://anaconda.org/maxibor/adrsm) [![Build Status](https://travis-ci.org/maxibor/adrsm.svg?branch=master)](https://travis-ci.org/maxibor/adrsm)

* * *

<img src="./img/logo_adrsm.png" width="300">

# Introduction

ADRSM (Ancient DNA Read Simulator for Metagenomics) is a tool designed to simulate the paired-end sequencing of a metagenomic community. ADRSM allows you to control precisely the amount of DNA from each organism in the community, which can be used to benchmark different metagenomics methods.

# Dependencies

-   [Conda](https://conda.io/miniconda.html)  

# Installation

    conda install -c maxibor adrsm

# Usage

    adrsm -d ./data/genomes ./data/short_genome_list.csv

# Output

-   `metagenome.{1,2}.fastq` : Simulated paired end reads
-   `stats.csv` : Statistics of simulated metagenome (organism, percentage of organism's DNA in metagenome)

# Help

    maxime@gph:~$ adrsm --help
    usage: ADRSM [-h] [-d DIRECTORY] [-r READLENGTH] [-l LENSTDEV] [-fwd FWDADAPT]
                 [-rev REVADAPT] [-e ERROR] [-p GEOM_P] [-m MIN] [-M MAX]
                 [-o OUTPUT] [-q QUALITY] [-s STATS]
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
      -p GEOM_P      Geometric distribution parameter for deamination. Default =
                     0.5
      -m MIN         Deamination substitution base frequency. Default = 0.001
      -M MAX         Deamination substitution max frequency. Default = 0.3
      -o OUTPUT      Output file basename. Default = ./metagenome.*
      -q QUALITY     Base quality encoding. Default = d (PHRED+64)
      -s STATS       Statistic file. Default = stats.csv

## Genome directory

Each genome `fasta` file must be named after the name of the organism. (example: [data/genomes](./data/genomes))

## Configuration file (`confFile`)

The configuration file is a `.csv` file describing, one line per genome, the mean insert size, and the expected genome coverage.
Example [short_genome_list.csv](./data/short_genome_list.csv):

| genome (mandatory)                  | insert_size (mandatory) | coverage (mandatory) | deamination (mandatory) |
| ----------------------------------- | ----------------------- | -------------------- | ----------------------- |
| Agrobacterium_tumefaciens_genome.fa | 47                      | 0.1                  | yes                     |
| Bacillus_anthracis_genome.fa        | 48                      | 0.2                  | no                      |

## Note on Coverage

Given the sequencing error, and the random choice of inserts, the target coverage might differ slightly from the real coverage (fig 1)  

<img src="./img/coverage_plot.png" width="600">  

**Figure 1:** Coverage plot for simulated sequencing of _Elephas maximus_ mitocondria. Aligned with Bowtie2 (default-parameters). Read-length = 76, insert-length = 200.

## Note on Deamination simulation

The deamination is modeled using a [Geometric distribution](https://en.wikipedia.org/wiki/Geometric_distribution)
With the default parameters, the substitution frequency is depicted in fig 2:  

<img src="./img/geometric_model.png" width="300">  

**Figure 2:** Substitution frequency.

For each nucleotide, a random number `Pu` is sampled from an <a href="https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)">uniform distribution</a> (of support [0 ,1]) and compared to the corresponding value `Pg` of the rescaled geometric PMF at this nucleotide.  
If `Pg >= Pu`, the base is substituted (fig 3).

<img src="./img/geometric_distribution.png" width="300">  

**Figure 3:** Substitutions distribution along a DNA insert, with default parameters.

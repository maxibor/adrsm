#!/usr/bin/env python

from numpy import random as npr
from .lib import adrsmlib as ad
from . import __version__
import click

@click.command()
@click.version_option(__version__)
@click.argument('confFile', type=click.Path(exists=True, 
                                            readable=True, 
                                            resolve_path=True))
@click.option('-r',
              '--readLength',
              default='76',
              type=int,
              show_default=True,
              help='Average read length')
@click.option('-n',
              '--nbinom',
              default=8,
              type=int,
              show_default=True,
              help='n parameter for Negative Binomial insert length distribution')
@click.option('-fwd',
              '--fwdAdapt',
              default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG',
              type=str,
              show_default=True,
              help='Forward adaptor sequence')
@click.option('-rev',
              '--revAdapt',
              default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
              type=str,
              show_default=True,
              help='Reverse adaptor sequence')
@click.option('-p', 
              '--geom_p',
              default=0.5,
              type=click.FloatRange(min=0.0, max=1.0),
              show_default=True,
              help='Geometric distribution parameter for deamination')
@click.option('-m',
              '--minD',
              default=0.01,
              type=click.FloatRange(min=0.0, max=1.0),
              show_default=True,
              help='Deamination substitution base frequency')   
@click.option('-M',
              '--maxD',
              default=0.3,
              type=click.FloatRange(min=0.0, max=1.0),
              show_default=True,
              help='Deamination substitution max frequency')
@click.option('-s',
              '--seed',
              default=42,
              type=int,
              show_default=True,
              help='Seed for random generator generator')
@click.option('-t',
              '--threads',
              default=2,
              type=click.IntRange(min=1, max=1024),
              show_default=True,
              help='Number of threads for parallel processing')
@click.option('-o',
              '--output',
              default='./metagenome',
              type=click.Path(file_okay=True, writable=True, resolve_path=True),
              show_default=True,
              help='Fastq output file basename')
@click.option('-s',
              '--stats',
              default='./stats.csv',
              type=click.Path(file_okay=True, writable=True, resolve_path=True),
              show_default=True,
              help='Summary statistics file')


def cli(no_args_is_help=True, **kwargs):
    """\b
    ==================================================
    ADRSM: Ancient DNA Read Simulator for Metagenomics
    Author: Maxime Borry
    Contact: <borry[at]shh.mpg.de>
    Homepage & Documentation: github.com/maxibor/adrsm

    CONFFILE: path to ADRSM configuration file
    """
    main(**kwargs)

def read_config(infile):
    """
    READS CONFIG FILE AND RETURNS CONFIG DICT
    """
    genomes = {}
    with open(infile, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip()
            splitline = line.split(",")
            agenome = splitline[0].replace(" ", "")
            ainsert = int(splitline[1].replace(" ", ""))
            acov = float(splitline[2].replace(" ", ""))
            deambool = str(splitline[3].replace(" ", ""))
            deamination = ad.parse_yes_no(deambool)
            if len(splitline) > 4 and float(splitline[4].replace(" ", "")) != 0.0:
                mutate = True
                mutrate = float(splitline[4].replace(" ", ""))
                age = float(splitline[5].replace(" ", ""))
            else:
                mutate = False
                mutrate = 0
                age = 0

            genomes[agenome] = {'size': ainsert,
                                'cov': acov, 'deam': deamination, 'mutate': mutate, 'mutrate': mutrate, 'age': age}
    return(genomes)


def main(conffile, readlength, nbinom, fwdadapt, revadapt, geom_p, mind, maxd, seed, threads, output, stats):
    MINLENGTH = 20
    npr.seed(seed)
    fastq_list = []
    stat_dict = {}
    all_genomes = read_config(conffile)
    for agenome in all_genomes.keys():
        stat_and_run = ad.run_read_simulation_multi(INFILE=agenome,
                                                    COV=all_genomes[agenome]['cov'],
                                                    READLEN=readlength,
                                                    INSERLEN=all_genomes[agenome]['size'],
                                                    NBINOM=nbinom,
                                                    A1=fwdadapt,
                                                    A2=revadapt,
                                                    MINLENGTH=MINLENGTH,
                                                    MUTATE=all_genomes[agenome]['mutate'],
                                                    MUTRATE=all_genomes[agenome]['mutrate'],
                                                    AGE=all_genomes[agenome]['age'],
                                                    DAMAGE=all_genomes[agenome]['deam'],
                                                    GEOM_P=geom_p,
                                                    THEMIN=mind,
                                                    THEMAX=maxd,
                                                    PROCESS=threads,
                                                    FASTQ_OUT=output)
        stat_dict[ad.get_basename(agenome)] = stat_and_run

    ad.write_stat(stat_dict=stat_dict, stat_out=stats)
    print("\n-- ADRSM v" + __version__ +
          " finished generating this mock metagenome --")
    print("-- FASTQ files written to " + output +
          ".1.fastq and " + output + ".2.fastq --")
    print("-- Statistic file written to " + stats + " --")

if __name__ == "__main__":
    cli()
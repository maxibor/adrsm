#!/usr/bin/env python

from numpy import random as npr
from adrsm.lib import adrsmlib as ad
from adrsm import __version__
import click
import csv
import random


@click.command()
@click.version_option(__version__)
@click.argument(
    "confFile", type=click.Path(exists=True, readable=True, resolve_path=True)
)
@click.option(
    "-r",
    "--readLength",
    default="76",
    type=int,
    show_default=True,
    help="Average read length",
)
@click.option(
    "-n",
    "--nbinom",
    default=8,
    type=int,
    show_default=True,
    help="n parameter for Negative Binomial insert length distribution",
)
@click.option(
    "-fwd",
    "--fwdAdapt",
    default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    type=str,
    show_default=True,
    help="Forward adaptor sequence",
)
@click.option(
    "-rev",
    "--revAdapt",
    default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    type=str,
    show_default=True,
    help="Reverse adaptor sequence",
)
@click.option(
    "-p",
    "--geom_p",
    default=0.5,
    type=click.FloatRange(min=0.0, max=1.0),
    show_default=True,
    help="Geometric distribution parameter for deamination",
)
@click.option(
    "-m",
    "--minD",
    default=0.01,
    type=click.FloatRange(min=0.0, max=1.0),
    show_default=True,
    help="Deamination substitution base frequency",
)
@click.option(
    "-M",
    "--maxD",
    default=0.3,
    type=click.FloatRange(min=0.0, max=1.0),
    show_default=True,
    help="Deamination substitution max frequency",
)
@click.option(
    "-e",
    "--effort",
    default=100,
    type=int,
    show_default=True,
    help="Sequencing effort, maximum number of reads to be generated",
)
@click.option(
    "-s",
    "--seed",
    default=42,
    type=int,
    show_default=True,
    help="Seed for random generator generator",
)
@click.option(
    "-t",
    "--threads",
    default=2,
    type=click.IntRange(min=1, max=1024),
    show_default=True,
    help="Number of threads for parallel processing",
)
@click.option(
    "-o",
    "--output",
    default="./metagenome",
    type=click.Path(file_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="Fastq output file basename",
)
@click.option(
    "-s",
    "--stats",
    default="./stats.csv",
    type=click.Path(file_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="Summary statistics file",
)
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
        csv_file = csv.DictReader(f)
        for row in csv_file:
            row = dict(row)
            row_name = row["genome(mandatory)"]
            row.pop("genome(mandatory)")
            genomes[row_name] = row

    for genome in genomes:
        genomes[genome] = {k.replace(" ", ""): v for k, v in genomes[genome].items()}
        genomes[genome]["size"] = int(genomes[genome].pop("insert_size(mandatory)"))
        genomes[genome]["cov"] = float(genomes[genome].pop("coverage(mandatory)"))
        genomes[genome]["deam"] = (
            genomes[genome].pop("deamination(mandatory)").replace(" ", "")
        )

        genomes[genome]["mutrate"] = genomes[genome].pop("mutation_rate(optional)")
        if genomes[genome]["mutrate"] is None or "NA" in genomes[genome]["mutrate"]:
            genomes[genome]["mutate"] = False
            genomes[genome]["mutrate"] = 0
            genomes[genome]["age"] = 0
        else:
            genomes[genome]["mutate"] = True
            genomes[genome]["mutrate"] = float(genomes[genome]["mutrate"])
            genomes[genome]["age"] = int(genomes[genome].pop("age(optional)"))
    return genomes


def main(
    conffile,
    readlength,
    nbinom,
    fwdadapt,
    revadapt,
    geom_p,
    mind,
    maxd,
    effort,
    seed,
    threads,
    output,
    stats,
):
    MINLENGTH = 20
    npr.seed(seed)
    fastq_list = []
    stat_dict = {}
    all_genomes = read_config(conffile)
    all_reads = []
    print("-- ADRSM v" + __version__ + " --")
    for agenome in all_genomes.keys():
        reads, stat_and_run = ad.run_read_simulation_multi(
            INFILE=agenome,
            COV=all_genomes[agenome]["cov"],
            READLEN=readlength,
            INSERLEN=all_genomes[agenome]["size"],
            NBINOM=nbinom,
            A1=fwdadapt,
            A2=revadapt,
            MINLENGTH=MINLENGTH,
            MUTATE=all_genomes[agenome]["mutate"],
            MUTRATE=all_genomes[agenome]["mutrate"],
            AGE=all_genomes[agenome]["age"],
            DAMAGE=all_genomes[agenome]["deam"],
            GEOM_P=geom_p,
            THEMIN=mind,
            THEMAX=maxd,
            PROCESS=threads,
        )
        stat_dict[ad.get_basename(agenome)] = stat_and_run
        all_reads.extend(reads)

    if len(all_reads) > effort:
        all_reads = random.sample(all_reads, effort)
    ad.write_fastq_multi(all_reads, output)
    ad.write_stat(stat_dict=stat_dict, stat_out=stats)
    print(
        "\n-- ADRSM v"
        + __version__
        + " finished generating "
        + str(len(all_reads))
        + " reads for this mock metagenome --"
    )
    print(
        "-- FASTQ files written to "
        + output
        + ".1.fastq.gz and "
        + output
        + ".2.fastq.gz --"
    )
    print("-- Statistic file written to " + stats + " --")


if __name__ == "__main__":
    cli()
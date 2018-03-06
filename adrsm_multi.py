#!/usr/bin/env python

from numpy import random as npr
import lib.adrsmlib as ad

import argparse

def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
    prog='MetaBenReadSim',
    description='Metagenomic Benchmarking Read Simulator for ancient DNA')
    parser.add_argument('confFile', help="path to configuration file")
    parser.add_argument(
    '-d',
    dest="directory",
    default = ".",
    help="path to genome directory. Default = .")
    parser.add_argument(
    '-r',
    dest='readLength',
    default=76,
    help="Average read length. Default = 76")
    parser.add_argument(
    '-l',
    dest="lenStdev",
    default=10,
    help="Insert length standard deviation. Default = 10")
    parser.add_argument(
    '-fwd',
    dest="fwdAdapt",
    default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    help="Forward adaptor. Default = AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
    parser.add_argument(
    '-rev',
    dest="revAdapt",
    default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    help="Reverse adaptor. Default = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
    parser.add_argument(
    '-e',
    dest="error",
    default=0.01,
    help="Illumina sequecing error. Default = 0.01")
    parser.add_argument(
    '-o',
    dest="output",
    default="metagenome",
    help="Output file basename. Default = ./metagenome.*")
    parser.add_argument(
    '-s',
    dest="stats",
    default="stats.csv",
    help="Statistic file. Default = stats.csv")

    args = parser.parse_args()

    infile = args.confFile
    gendir = args.directory
    readlen = int(args.readLength)
    lendev = int(args.lenStdev)
    a1 = args.fwdAdapt
    a2 = args.revAdapt
    err = float(args.error)
    outfile= args.output
    stats = args.stats

    return(infile, gendir, readlen, lendev, a1, a2, err, outfile, stats)

def read_config(infile, gendir):
    """
    READS CONFIG FILE AND RETURNS CONFIG DICT
    """
    genomes = {}
    with open(infile, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip()
            splitline = line.split(",")
            agenome = splitline[0].replace(" ","")
            ainsert = int(splitline[1].replace(" ",""))
            acov = float(splitline[2].replace(" ",""))
            genomes[gendir+"/"+agenome] = [ainsert, acov]
    return(genomes)

if __name__ == "__main__":
    INFILE, GENDIR, READLEN, LENDEV, A1, A2, ERR, OUTFILE, STATS = _get_args()
    MINLENGTH = 20

    genome_dict = {}
    stat_dict = {}
    all_genomes = read_config(INFILE, GENDIR)
    for agenome in all_genomes.keys():
        stat_and_run = ad.run_read_simulation_multi(INFILE = agenome,
                                  NREAD = None,
                                  COV = all_genomes[agenome][1],
                                  READLEN = READLEN,
                                  INSERLEN = all_genomes[agenome][0] ,
                                  LENDEV = LENDEV,
                                  A1 = A1,
                                  A2 = A2,
                                  MINLENGTH = MINLENGTH,
                                  ERR = ERR,
                                  fastq_dict = genome_dict)
        stat_dict[ad.get_basename(agenome)] = stat_and_run

    ad.write_fastq_multi(fastq_dict=genome_dict, outputfile=OUTFILE)
    ad.write_stat(stat_dict=stat_dict, stat_out=STATS)

#!/usr/bin/env python

import sys
import requests
import os
import subprocess
from numpy import random as npr
import multiprocessing
import pickle
from functools import partial
from pkg_resources import resource_filename
from . import sequencefunctions as sf
from . import markov as mk
from xopen import xopen


def parse_yes_no(astring):
    if "yes" in astring:
        return True
    elif "no" in astring:
        return False
    else:
        sys.exit("Please specify deamination (yes | no)")


def get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return basename


def read_fasta(file_name):
    """
    READS FASTA FILE, RETURNS SEQUENCE AS STRING
    INPUT:
        file_name(string): path to fasta file
    OUPUT:
        result(string): all of the sequences in fasta file, concatenated
    """
    result = ""
    # fastadict = {}
    if file_name.endswith(".gz"):
        with xopen(file_name, "r") as f:
            for line in f:
                if line[0] == ">":
                    # seqname = line[1:]
                    # fastadict[seqname] = []
                    continue
                else:
                    line = line.rstrip()
                    # fastadict[seqname].append(line)
                    result += line
    else:
        with open(file_name, "r") as f:
            for line in f:
                if line[0] == ">":
                    # seqname = line[1:]
                    # fastadict[seqname] = []
                    continue
                else:
                    line = line.rstrip()
                    # fastadict[seqname].append(line)
                    result += line
    return [result, len(result)]


def write_fastq_multi(fastq_list, outputfile, compressed=True):
    if compressed:
        with xopen(outputfile + ".1.fastq.gz", "ab") as f1:
            with xopen(outputfile + ".2.fastq.gz", "ab") as f2:
                for read in fastq_list:
                    f1.write(read[0].encode())
                    f2.write(read[1].encode())
    else:
        with open(outputfile + ".1.fastq", "a") as f1:
            with open(outputfile + ".2.fastq", "a") as f2:
                for read in fastq_list:
                    f1.write(read[0])
                    f2.write(read[1])


def markov_wrapper_fwd(times):
    a = 0
    while a == 0:
        a = mk.mchain(
            starts=MARKOV_START_FWD,
            kmers=MARKOV_DICT_FWD,
            readsize=READSIZE,
            order=MARKOV_ORDER,
        )
    return a


def markov_wrapper_rev(times):
    a = 0
    while a == 0:
        a = mk.mchain(
            starts=MARKOV_START_REV,
            kmers=MARKOV_DICT_REV,
            readsize=READSIZE,
            order=MARKOV_ORDER,
        )
    return a


def markov_multi_fwd(process, nreads):
    myIter = range(nreads)
    with multiprocessing.Pool(process) as p:
        r = p.map(markov_wrapper_fwd, myIter)
    return r


def markov_multi_rev(process, nreads):
    myIter = range(nreads)
    with multiprocessing.Pool(process) as p:
        r = p.map(markov_wrapper_rev, myIter)
    return r


def get_fwd_qual():
    try:
        ret = pickle.load(open("data/quality/fwd_qual.p", "rb"))
        return ret
    except FileNotFoundError:
        path = resource_filename("adrsm", "/data/quality/fwd_qual.p")
        ret = pickle.load(open(path, "rb"))
        return ret


def get_rev_qual():
    try:
        ret = pickle.load(open("data/quality/fwd_qual.p", "rb"))
        return ret
    except FileNotFoundError:
        path = resource_filename("adrsm", "/data/quality/rev_qual.p")
        ret = pickle.load(open(path, "rb"))
        return ret


def multi_run(
    iterables,
    name,
    mutate,
    mutrate,
    damage,
    geom_p,
    themin,
    themax,
    fwd_adaptor,
    rev_adaptor,
    read_length,
    process,
):
    partial_run = partial(
        sf.generate_fq,
        name=name,
        mutate=mutate,
        mutrate=mutrate,
        damage=damage,
        geom_p=geom_p,
        themin=themin,
        themax=themax,
        fwd_adaptor=fwd_adaptor,
        rev_adaptor=rev_adaptor,
        read_length=read_length,
    )
    with multiprocessing.Pool(process) as p:
        r = p.map(partial_run, iterables)
    return r


def run_read_simulation_multi(
    INFILE,
    COV,
    READLEN,
    INSERLEN,
    NBINOM,
    A1,
    A2,
    MINLENGTH,
    MUTATE,
    MUTRATE,
    AGE,
    DAMAGE,
    GEOM_P,
    THEMIN,
    THEMAX,
    PROCESS,
):
    print("===================\n===================")
    print("Genome: ", INFILE)
    print("Coverage: ", COV)
    print("Read length: ", READLEN)
    print("Mean Insert length: ", INSERLEN)
    print("n parameter for Negative Binomial insert length distribution: ", NBINOM)
    print("Adaptor 1: ", A1)
    print("Adaptor 2: ", A2)
    print("Mutation rate (bp/year):", MUTRATE)
    print("Age (years):", AGE)
    print("Deamination:", DAMAGE)
    nread = None
    global READSIZE
    global MARKOV_ORDER
    global QUALIT_FWD
    global MARKOV_SEED_FWD
    global MARKOV_START_FWD
    global MARKOV_DICT_FWD
    global QUALIT_REV
    global MARKOV_SEED_REV
    global MARKOV_START_REV
    global MARKOV_DICT_REV

    READSIZE = READLEN

    basename = get_basename(INFILE)
    fasta = read_fasta(INFILE)

    nread = int((fasta[1] / INSERLEN) * COV)
    print("Number of reads: ", nread)
    print("-------------------")

    MARKOV_ORDER = 10
    QUALIT_FWD = get_fwd_qual()
    QUALIT_REV = get_rev_qual()
    MARKOV_SEED_FWD = mk.generate_kmer(
        qualities=QUALIT_FWD, order=MARKOV_ORDER, readsize=READLEN
    )
    MARKOV_SEED_REV = mk.generate_kmer(
        qualities=QUALIT_REV, order=MARKOV_ORDER, readsize=READLEN
    )
    MARKOV_START_FWD = MARKOV_SEED_FWD[0]
    MARKOV_START_REV = MARKOV_SEED_REV[0]
    MARKOV_DICT_FWD = MARKOV_SEED_FWD[1]
    MARKOV_DICT_REV = MARKOV_SEED_REV[1]

    # negative_binomial parameters
    prob = NBINOM / (NBINOM + INSERLEN)
    fragment_lengths = npr.negative_binomial(NBINOM, prob, nread)

    # Define Mutation rate
    if MUTATE:
        correct_mutrate = (MUTRATE * AGE) / fasta[1]
    else:
        correct_mutrate = 0

    # Prepare fragments and errors
    all_fragments = sf.random_insert(fasta, fragment_lengths, READLEN, MINLENGTH)
    fwd_illu_err = markov_multi_fwd(process=PROCESS, nreads=len(all_fragments))
    rev_illu_err = markov_multi_rev(process=PROCESS, nreads=len(all_fragments))

    runlist = sf.prepare_run(
        all_frag=all_fragments, all_fwd_err=fwd_illu_err, all_rev_err=rev_illu_err
    )

    result = multi_run(
        iterables=runlist,
        name=basename,
        mutate=MUTATE,
        mutrate=correct_mutrate,
        damage=DAMAGE,
        geom_p=GEOM_P,
        themin=THEMIN,
        themax=THEMAX,
        fwd_adaptor=A1,
        rev_adaptor=A2,
        read_length=READLEN,
        process=PROCESS,
    )

    # write_fastq_multi(fastq_list=result, outputfile=FASTQ_OUT)
    return result, [nread * INSERLEN, INSERLEN, COV, DAMAGE]


def specie_to_taxid(specie):
    """
    Takes a specie_name (ex: Mus_musculus), makes a call to JGI
    taxonomy API, and returns taxonomy id.

    INPUT:
        specie(string) ex: "Mus musculus"
    OUPUT:
        taxid(str) "10090"
    """

    request = "http://taxonomy.jgi-psf.org/tax/pt_name/" + specie
    response = requests.get(request)
    answer = response.text
    return answer


def write_stat(stat_dict, stat_out):
    nbases = []
    for akey in stat_dict:
        nbases.append(stat_dict[akey][0])
    totbases = sum(nbases)
    with open(stat_out, "w") as fs:
        fs.write(
            "Organism,taxonomy_id,percentage of metagenome,mean_insert_length,target_coverage,deamination\n"
        )
        for akey in stat_dict:
            taxid = specie_to_taxid(akey)
            fs.write(
                akey
                + ","
                + str(taxid)
                + ","
                + str(round(stat_dict[akey][0] / totbases, 2))
                + ","
                + str(stat_dict[akey][1])
                + ","
                + str(stat_dict[akey][2])
                + ","
                + str(stat_dict[akey][3])
                + "\n"
            )

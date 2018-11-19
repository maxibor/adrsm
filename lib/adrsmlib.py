#!/usr/bin/env python

import sys
import requests
import os
import subprocess
from numpy import random as npr
import multiprocessing
import pickle
from functools import partial
from . import sequencefunctions as sf
from . import markov as mk


def parse_yes_no(astring):
    if "yes" in astring:
        return(True)
    elif "no" in astring:
        return(False)
    else:
        sys.exit("Please specify deamination (yes | no)")


def get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)


def add_mutation_multi(sequences, mutrate, process):
    mutate_partial = partial(sf.mutate, mutrate=mutrate)
    print("Mutating...")
    with multiprocessing.Pool(process) as p:
        mutseq = p.map(mutate_partial, sequences)
    return(list(mutseq))


def reverse_complement_multi(all_inserts, process):
    print("Reverse Complemeting...")
    with multiprocessing.Pool(process) as p:
        r = p.map(sf.reverse_complement, all_inserts)
    return(r)


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
    with open(file_name, "r") as f:
        for line in f:
            if line[0] == ">":
                # seqname = line[1:]
                # fastadict[seqname] = []
                continue
            else:
                line = line.rstrip()
                # fastadict[seqname].append(line)
                result = result + line
    return([result, len(result)])


def complement_read_single(all_inserts, adaptor, read_length):
    result = []
    for insert in all_inserts:
        result.append(sf.complement_read(
            insert, adaptor=adaptor, read_length=read_length))
    return(result)


def complement_read_multi(all_inserts, adaptor, read_length, process):
    complement_read_partial = partial(
        sf.complement_read, adaptor=adaptor, read_length=read_length)
    with multiprocessing.Pool(process) as p:
        r = p.map(complement_read_partial, all_inserts)
    return(r)


def add_damage_single(all_inserts, geom_p, scale_min, scale_max):
    for i in range(0, len(all_inserts)):
        all_inserts[i] = sf.add_damage(
            insert=geom_p, scale_min=scale_min, scale_max=scale_max)
    return(all_inserts)


def add_damage_multi(all_inserts, geom_p, scale_min, scale_max, process):
    add_damage_partial = partial(
        sf.add_damage, geom_p=geom_p, scale_min=scale_min, scale_max=scale_max)
    print("Adding damage...")
    with multiprocessing.Pool(process) as p:
        r = p.map(add_damage_partial, all_inserts)
    return(r)


def add_error_single(all_reads, error_rate):
    for i in range(0, len(all_reads)):
        all_reads[i] = sf.add_error(read=i, error_rate=error_rate)
    return(all_reads)


def add_error_multi(all_reads, error_rate, process):
    add_error_partial = partial(sf.add_error, error_rate=error_rate)
    with multiprocessing.Pool(process) as p:
        r = p.map(add_error_partial, all_reads)
    return(r)


def prepare_fastq(fastq_dict, fwd_reads, rev_reads, fwd_illu_err, rev_illu_err, basename, read_length):
    fastq_dict[basename] = [[] for i in range(2)]
    cnt = 1
    for read1, read2, err1, err2 in zip(fwd_reads, rev_reads, fwd_illu_err, rev_illu_err):
        read1 = read1.rstrip()
        read2 = read2.rstrip()
        towrite_fwd = "@" + basename + "_" + \
            str(cnt) + "/1" + "\n" + read1 + \
            "\n+\n" + err1 + "\n"
        fastq_dict[basename][0].append(towrite_fwd)
        towrite_rev = "@" + basename + "_" + \
            str(cnt) + "/2" + "\n" + read2 + \
            "\n+\n" + err2 + "\n"
        fastq_dict[basename][1].append(towrite_rev)
        cnt += 1
    return(fastq_dict)


def write_fastq_multi(fastq_dict, outputfile):
    with open(outputfile + ".1.fastq", "w") as f1:
        with open(outputfile + ".2.fastq", "w") as f2:
            for akey in fastq_dict.keys():
                for reads1 in fastq_dict[akey][0]:
                    f1.write(reads1)
                for reads2 in fastq_dict[akey][1]:
                    f2.write(reads2)


def markov_wrapper_fwd(times):
    a = 0
    while(a == 0):
        a = mk.mchain(starts=MARKOV_START_FWD, kmers=MARKOV_DICT_FWD,
                      readsize=READSIZE, order=MARKOV_ORDER)
    return(a)


def markov_wrapper_rev(times):
    a = 0
    while(a == 0):
        a = mk.mchain(starts=MARKOV_START_REV, kmers=MARKOV_DICT_REV,
                      readsize=READSIZE, order=MARKOV_ORDER)
    return(a)


def markov_multi_fwd(process, nreads):
    myIter = range(nreads)
    with multiprocessing.Pool(process) as p:
        r = p.map(markov_wrapper_fwd, myIter)
    return(r)


def markov_multi_rev(process, nreads):
    myIter = range(nreads)
    with multiprocessing.Pool(process) as p:
        r = p.map(markov_wrapper_rev, myIter)
    return(r)


def get_fwd_qual():
    try:
        ret = pickle.load(open("data/quality/fwd_qual.p", 'rb'))
        return(ret)
    except FileNotFoundError:
        cmd = "which adrsm"
        res = subprocess.check_output(cmd, shell=True)
        res = res.decode('utf-8').rstrip()
        path = "/".join(res.split("/")[:-2])+"/data/quality/fwd_qual.p"
        ret = pickle.load(open(path, 'rb'))
        return(ret)


def get_rev_qual():
    try:
        ret = pickle.load(open("data/quality/fwd_qual.p", 'rb'))
        return(ret)
    except FileNotFoundError:
        cmd = "which adrsm"
        res = subprocess.check_output(cmd, shell=True)
        res = res.decode('utf-8').rstrip()
        path = "/".join(res.split("/")[:-2])+"/data/quality/rev_qual.p"
        ret = pickle.load(open(path, 'rb'))
        return(ret)


def run_read_simulation_multi(INFILE, COV, READLEN, INSERLEN, NBINOM, A1, A2, MINLENGTH, MUTATE, MUTRATE, AGE, ERR,  DAMAGE, GEOM_P, THEMIN, THEMAX, fastq_dict, PROCESS):
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
    print("Sequencing Error rate", ERR)
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

    MARKOV_ORDER = 7
    QUALIT_FWD = get_fwd_qual()
    QUALIT_REV = get_rev_qual()
    MARKOV_SEED_FWD = mk.generate_kmer(
        qualities=QUALIT_FWD, order=MARKOV_ORDER, readsize=READLEN)
    MARKOV_SEED_REV = mk.generate_kmer(
        qualities=QUALIT_REV, order=MARKOV_ORDER, readsize=READLEN)
    MARKOV_START_FWD = MARKOV_SEED_FWD[0]
    MARKOV_START_REV = MARKOV_SEED_REV[0]
    MARKOV_DICT_FWD = MARKOV_SEED_FWD[1]
    MARKOV_DICT_REV = MARKOV_SEED_REV[1]

    # negative_binomial parameters
    prob = NBINOM / (NBINOM + INSERLEN)
    insert_lengths = npr.negative_binomial(NBINOM, prob, nread)

    all_inserts = sf.random_insert(fasta, insert_lengths, READLEN, MINLENGTH)

    if MUTATE:
        correct_mutrate = (MUTRATE * AGE) / fasta[1]
        all_inserts = add_mutation_multi(
            sequences=all_inserts, mutrate=correct_mutrate, process=PROCESS)

    if DAMAGE:
        all_inserts = add_damage_multi(
            all_inserts=all_inserts,
            geom_p=GEOM_P,
            scale_min=THEMIN,
            scale_max=THEMAX,
            process=PROCESS)

    fwd_inserts = all_inserts
    # rev_inserts = [sf.reverse_complement(i) for i in all_inserts]
    rev_inserts = reverse_complement_multi(
        all_inserts=all_inserts, process=PROCESS)
    print("Adding adaptors to forward read...")
    # fwd_reads = complement_read_single(fwd_inserts, A1, READLEN)
    fwd_reads = complement_read_multi(fwd_inserts, A1, READLEN, PROCESS)
    print("Adding sequencing error to forward read...")
    # fwd_reads = add_error_single(fwd_reads, ERR)
    fwd_reads = add_error_multi(fwd_reads, ERR, PROCESS)
    fwd_illu_err = markov_multi_fwd(process=PROCESS, nreads=len(fwd_reads))
    print("Adding adaptors to reverse read...")
    # rev_reads = complement_read_single(rev_inserts, A2, READLEN)
    rev_reads = complement_read_multi(rev_inserts, A2, READLEN, PROCESS)
    print("Adding sequencing error to reverse read...")
    # rev_reads = add_error_single(rev_reads, ERR)
    rev_reads = add_error_multi(rev_reads, ERR, PROCESS)
    rev_illu_err = markov_multi_rev(process=PROCESS, nreads=len(rev_reads))

    prepare_fastq(fastq_dict=fastq_dict,
                  fwd_reads=fwd_reads,
                  rev_reads=rev_reads,
                  fwd_illu_err=fwd_illu_err,
                  rev_illu_err=rev_illu_err,
                  basename=basename,
                  read_length=READLEN
                  )
    return([nread * INSERLEN, INSERLEN, COV, DAMAGE])


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
    return(answer)


def write_stat(stat_dict, stat_out):
    nbases = []
    for akey in stat_dict:
        nbases.append(stat_dict[akey][0])
    totbases = sum(nbases)
    with open(stat_out, "w") as fs:
        fs.write(
            "Organism,taxonomy_id,percentage of metagenome,mean_insert_length,target_coverage,deamination\n")
        for akey in stat_dict:
            taxid = specie_to_taxid(akey)
            fs.write(akey + "," + str(taxid) + "," + str(round(stat_dict[akey][0] / totbases, 2)) + "," + str(
                stat_dict[akey][1]) + "," + str(stat_dict[akey][2]) + "," + str(stat_dict[akey][3]) + "\n")

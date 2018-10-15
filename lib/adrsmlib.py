#!/usr/bin/env python

import sys
import requests
import numpy as np
from numpy import random as npr
import multiprocessing
from functools import partial
from scipy.stats import geom


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


def scale(x, themin, themax):
    return(np.interp(x, (x.min(), x.max()), (themin, themax)))


def reverse_complement(dna):
    dna = dna.upper()
    '''
    Reverse complement a DNA string
    '''
    dna = dna[::-1]
    revcom = []
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    for letter in dna:
        for key in complement.keys():
            if letter == key:
                revcom.append(complement[key])

    return "".join(revcom)


def reverse_complement_multi(all_inserts, process):
    with multiprocessing.Pool(process) as p:
        r = list(p.map(reverse_complement, all_inserts))
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
    with open(file_name, "r") as f:
        for line in f:
            if not line.startswith(">"):
                line = line.rstrip()
                result = result + line
    return([result, len(result)])


def random_insert(read_fasta_out, insert_lengths, read_length, minlen):
    genome = read_fasta_out[0]
    genome_length = read_fasta_out[1]
    result = []
    for i in insert_lengths:
        if i >= minlen:
            insert_start = npr.randint(0, genome_length - read_length)
            insert_end = insert_start + i + 1
            insert = genome[insert_start:insert_end]
            result.append(insert)
    return(result)


def _complement_read(insert, adaptor, read_length):
    inlen = len(insert)
    if inlen < read_length:
        diff = read_length - inlen
        to_add = adaptor[0:diff]
        read = insert + to_add
    elif inlen == read_length:
        read = insert
    elif inlen > read_length:
        read = insert[0:read_length]
    if len(read) == read_length:
        read = read.upper()
        read = list(read)
        for j in range(0, len(read)):
            if read[j] not in ["A", "T", "G", "C", "N"]:
                read[j] = "N"
        return("".join(read))


def complement_read_multi(all_inserts, adaptor, read_length, process):
    complement_read_partial = partial(
        _complement_read, adaptor=adaptor, read_length=read_length)
    with multiprocessing.Pool(process) as p:
        r = list(p.map(complement_read_partial, all_inserts))
    return(r)


def complement_read(all_inserts, adaptor, read_length):
    result = []
    for insert in all_inserts:
        inlen = len(insert)
        if inlen < read_length:
            diff = read_length - inlen
            to_add = adaptor[0:diff]
            read = insert + to_add
        elif inlen == read_length:
            read = insert
        elif inlen > read_length:
            read = insert[0:read_length]
        if len(read) == read_length:
            read = read.upper()
            read = list(read)
            for j in range(0, len(read)):
                if read[j] not in ["A", "T", "G", "C", "N"]:
                    read[j] = "N"
            result.append("".join(read))
    return(result)


def _add_damage(insert, geom_p, scale_min, scale_max):
    insert = list(insert)
    insertlen = len(insert)
    x = np.arange(1, insertlen + 1)
    geom_dist = scale(geom.pmf(x, geom_p), scale_min, scale_max)

    for j in range(0, insertlen):
        pos = j
        opp_pos = insertlen - 1 - j

        # C -> T deamination
        if insert[pos] == "C" and geom_dist[j] >= npr.rand():
            insert[pos] = "T"

        # G -> A deamination
        if insert[opp_pos] == "G" and geom_dist[j] >= npr.rand():
            insert[opp_pos] = "A"
    return("".join(insert))


def add_damage_multi(all_inserts, geom_p, scale_min, scale_max, process):
    add_damage_partial = partial(
        _add_damage, geom_p=geom_p, scale_min=scale_min, scale_max=scale_max)
    with multiprocessing.Pool(process) as p:
        r = list(p.map(add_damage_partial, all_inserts))
    return(r)


def add_damage(all_inserts, geom_p, scale_min, scale_max):
    for i in range(0, len(all_inserts)):
        insert = list(all_inserts[i])
        insertlen = len(insert)
        x = np.arange(1, insertlen + 1)
        geom_dist = scale(geom.pmf(x, geom_p), scale_min, scale_max)

        for j in range(0, insertlen):
            pos = j
            opp_pos = insertlen - 1 - j

            # C -> T deamination
            if insert[pos] == "C" and geom_dist[j] >= npr.rand():
                insert[pos] = "T"

            # G -> A deamination
            if insert[opp_pos] == "G" and geom_dist[j] >= npr.rand():
                insert[opp_pos] = "A"
        all_inserts[i] = "".join(insert)
    return(all_inserts)


def _add_error(read, error_rate):
    read = list(read)
    for j in range(0, len(read)):
        if read[j].upper() not in ["A", "T", "G", "C", "N"]:
            read[j] = "N"
        if npr.random() < error_rate:
            read[j] = npr.choice(["A", "T", "G", "C"])
    return("".join(read))


def add_error_multi(all_reads, error_rate, process):
    add_error_partial = partial(_add_error, error_rate=error_rate)
    with multiprocessing.Pool(process) as p:
        r = list(p.map(add_error_partial, all_reads))
    return(r)


def add_error(all_reads, error_rate):
    for i in range(0, len(all_reads)):
        read = list(all_reads[i])
        for j in range(0, len(read)):
            if read[j].upper() not in ["A", "T", "G", "C", "N"]:
                read[j] = "N"
            if npr.random() < error_rate:
                read[j] = npr.choice(["A", "T", "G", "C"])
        all_reads[i] = "".join(read)
    return(all_reads)


def prepare_fastq(fastq_dict, fwd_reads, rev_reads, basename, read_length, quality):
    fastq_dict[basename] = [[] for i in range(2)]
    cnt = 1
    for read1, read2 in zip(fwd_reads, rev_reads):
        read1 = read1.rstrip()
        read2 = read2.rstrip()
        readlen1 = len(read1)
        readlen2 = len(read2)
        towrite_fwd = "@" + basename + "_" + \
            str(cnt) + "/1" + "\n" + read1 + \
            "\n+\n" + quality * readlen1 + "\n"
        fastq_dict[basename][0].append(towrite_fwd)
        towrite_rev = "@" + basename + "_" + \
            str(cnt) + "/2" + "\n" + read2 + \
            "\n+\n" + quality * readlen2 + "\n"
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


def run_read_simulation_multi(INFILE, COV, READLEN, INSERLEN, NBINOM, A1, A2, MINLENGTH, ERR,  DAMAGE, GEOM_P, THEMIN, THEMAX, fastq_dict, QUALITY, PROCESS):
    print("===================")
    print("Genome: ", INFILE)
    print("Coverage: ", COV)
    print("Read length: ", READLEN)
    print("Mean Insert length: ", INSERLEN)
    print("n parameter for Negative Binomial insert length distribution: ", NBINOM)
    print("Adaptor 1: ", A1)
    print("Adaptor 2: ", A2)
    print("Quality :", QUALITY)
    print("Deamination:", DAMAGE)
    nread = None

    basename = get_basename(INFILE)
    fasta = read_fasta(INFILE)

    nread = int((fasta[1] / INSERLEN) * COV)
    print("Number of reads: ", nread)
    print("===================\n")

    # negative_binomial parameters
    prob = NBINOM / (NBINOM + INSERLEN)
    insert_lengths = npr.negative_binomial(NBINOM, prob, nread)

    all_inserts = random_insert(fasta, insert_lengths, READLEN, MINLENGTH)
    if DAMAGE:
        # all_inserts = add_damage(
        #     all_inserts=all_inserts,
        #     geom_p=GEOM_P,
        #     scale_min=THEMIN,
        #     scale_max=THEMAX)
        all_inserts = add_damage_multi(
            all_inserts=all_inserts,
            geom_p=GEOM_P,
            scale_min=THEMIN,
            scale_max=THEMAX,
            process=PROCESS)

    fwd_inserts = all_inserts
    # rev_inserts = [reverse_complement(i) for i in all_inserts]
    rev_inserts = reverse_complement_multi(
        all_inserts=all_inserts, process=PROCESS)
    # fwd_reads = complement_read(fwd_inserts, A1, READLEN)
    fwd_reads = complement_read_multi(fwd_inserts, A1, READLEN, PROCESS)
    # fwd_reads = add_error(fwd_reads, ERR)
    fwd_reads = add_error_multi(fwd_reads, ERR, PROCESS)
    # rev_reads = complement_read(rev_inserts, A2, READLEN)
    rev_reads = complement_read_multi(rev_inserts, A2, READLEN, PROCESS)
    # rev_reads = add_error(rev_reads, ERR)
    rev_reads = add_error_multi(rev_reads, ERR, PROCESS)

    prepare_fastq(fastq_dict=fastq_dict,
                  fwd_reads=fwd_reads,
                  rev_reads=rev_reads,
                  basename=basename,
                  read_length=READLEN,
                  quality=QUALITY)
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

import numpy as np
from numpy import random as npr
from scipy.stats import geom
import random


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


def complement_read(insert, adaptor, read_length):
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


def add_damage(insert, geom_p, scale_min, scale_max):
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


def add_error(read, error_rate):
    read = list(read)
    for j in range(0, len(read)):
        if read[j].upper() not in ["A", "T", "G", "C", "N"]:
            read[j] = "N"
        if npr.random() < error_rate:
            read[j] = npr.choice(["A", "T", "G", "C"])
    return("".join(read))


def mutate(nucleotide, mutrate):
    """
    alpha: Transitions
    beta: Transversions
    https://en.wikipedia.org/wiki/Mutation_rate
    """
    alpha = 0.2
    beta = 0.4
    a = int(10 * alpha)
    b = int(10 * beta)
    dmut = {'A': b * ['C'] + b * ['T'] + a * ['G'], 'C': b * ['A'] + b * ['G'] + a * [
        'T'], 'G': b * ['C'] + b * ['T'] + a * ['A'], 'T': b * ['A'] + b * ['G'] + a * ['C']}
    if npr.random() <= mutrate:
        new_nucl = random.choice(dmut[nucleotide])
        return(new_nucl)
    return(nucleotide)

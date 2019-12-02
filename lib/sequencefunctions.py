import numpy as np
from numpy import random as npr
from scipy.stats import geom
import random
from . import quality


def random_insert(read_fasta_out, insert_lengths, read_length, minlen):
    genome = read_fasta_out[0]
    genome_length = read_fasta_out[1]
    result = []
    for i in insert_lengths:
        if i >= minlen:
            insert_start = npr.randint(0, genome_length - i - 1)
            insert_end = insert_start + i + 1
            insert = genome[insert_start:insert_end]
            result.append(insert)
    return(result)


def scale(x, themin, themax):
    return(np.interp(x, (x.min(), x.max()), (themin, themax)))


class fragment ():

    def __init__(self, sequence, num, name, fwd_phred, rev_phred):
        """
        sequence(str) nucleotide sequence of fragment
        num(int) id of the sequence
        name(str) name of the genome
        """
        self.seq = sequence.upper()
        self.id = num
        self.name = name
        self.fwd_phred = fwd_phred
        self.rev_phred = rev_phred

    def __repr__(self):
        return(f"A fragment class object of sequence {self.seq}")

    def reverse_complement(self):
        '''
        Reverse complement a DNA string
        '''
        dna = self.seq[::-1]
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        revcom = [complement[i] for i in dna]
        self.revcom = "".join(revcom)
        return(self.revcom)

    def mutate_fwd(self, mutrate, alpha=0.4, beta=0.2):
        """
        alpha: Transitions
        beta: Transversions
        https://en.wikipedia.org/wiki/Mutation_rate
        """
        a = int(10 * alpha)
        b = int(10 * beta)
        newseq = ""
        dmut = {'A': b * ['C'] + b * ['T'] + a * ['G'], 'C': b * ['A'] + b * ['G'] + a * [
            'T'], 'G': b * ['C'] + b * ['T'] + a * ['A'], 'T': b * ['A'] + b * ['G'] + a * ['C']}
        for nuc in self.seq:
            if npr.random() <= mutrate:
                new_nucl = random.choice(dmut[nuc])
                newseq += new_nucl
            else:
                newseq += nuc
        return(self.seq)

    def mutate_rev(self, mutrate, alpha=0.4, beta=0.2):
        """
        alpha: Transitions
        beta: Transversions
        https://en.wikipedia.org/wiki/Mutation_rate
        """
        a = int(10 * alpha)
        b = int(10 * beta)
        newseq = ""
        dmut = {'A': b * ['C'] + b * ['T'] + a * ['G'], 'C': b * ['A'] + b * ['G'] + a * [
            'T'], 'G': b * ['C'] + b * ['T'] + a * ['A'], 'T': b * ['A'] + b * ['G'] + a * ['C']}
        for nuc in self.revcom:
            if npr.random() <= mutrate:
                new_nucl = random.choice(dmut[nuc])
                newseq += new_nucl
            else:
                newseq += nuc
        return(self.revcom)

    def add_damage_fwd(self, geom_p, scale_min, scale_max):

        insert = list(self.seq)
        insertlen = len(self.seq)
        x = np.arange(1, insertlen + 1)
        geom_dist = geom.pmf(x, geom_p)
        # geom_dist = scale(geom.pmf(x, geom_p), scale_min, scale_max)

        for j in range(0, insertlen):
            pos = j
            opp_pos = insertlen - 1 - j

            # C -> T deamination - deamination
            if npr.rand() <= scale_max:
                if insert[pos] == "C" and geom_dist[j] >= npr.rand():
                    insert[pos] = "T"

            # C -> T deamination - baseline
            if npr.rand() <= scale_min:
                if insert[pos] == "C":
                    insert[pos] = "T"

            # G -> A deamination
            if npr.rand() <= scale_max:
                if insert[opp_pos] == "G" and geom_dist[j] >= npr.rand():
                    insert[opp_pos] = "A"

            # C -> T deamination - baseline
            if npr.rand() <= scale_min:
                if insert[pos] == "G":
                    insert[pos] = "A"

        self.seq = "".join(insert)
        return(self.seq)

    def add_damage_rev(self, geom_p, scale_min, scale_max):

        insert = list(self.revcom)
        insertlen = len(self.revcom)
        x = np.arange(1, insertlen + 1)
        geom_dist = geom.pmf(x, geom_p)
        # geom_dist = scale(geom.pmf(x, geom_p), scale_min, scale_max)

        for j in range(0, insertlen):
            pos = j
            opp_pos = insertlen - 1 - j

           # C -> T deamination - deamination
            if npr.rand() <= scale_max:
                if insert[pos] == "C" and geom_dist[j] >= npr.rand():
                    insert[pos] = "T"

            # C -> T deamination - baseline
            if npr.rand() <= scale_min:
                if insert[pos] == "C":
                    insert[pos] = "T"

            # G -> A deamination
            if npr.rand() <= scale_max:
                if insert[opp_pos] == "G" and geom_dist[j] >= npr.rand():
                    insert[opp_pos] = "A"

            # C -> T deamination - baseline
            if npr.rand() <= scale_min:
                if insert[pos] == "G":
                    insert[pos] = "A"
                    
        self.revcom = "".join(insert)
        return(self.revcom)

    def complement_read_fwd(self, adaptor, read_length):
        inlen = len(self.seq)
        if inlen < read_length:
            diff = read_length - inlen
            to_add = adaptor[0:diff]
            read = self.seq + to_add
        elif inlen == read_length:
            read = self.seq
        elif inlen > read_length:
            read = self.seq[0:read_length]
        if len(read) == read_length:
            read = read.upper()
            read = list(read)
            for j in range(0, len(read)):
                if read[j] not in ["A", "T", "G", "C", "N"]:
                    read[j] = "N"
        self.fwd_read = "".join(read)
        return(self.fwd_read)

    def complement_read_rev(self, adaptor, read_length):
        inlen = len(self.revcom)
        if inlen < read_length:
            diff = read_length - inlen
            to_add = adaptor[0:diff]
            read = self.revcom + to_add
        elif inlen == read_length:
            read = self.revcom
        elif inlen > read_length:
            read = self.revcom[0:read_length]
        if len(read) == read_length:
            read = read.upper()
            read = list(read)
            for j in range(0, len(read)):
                if read[j] not in ["A", "T", "G", "C", "N"]:
                    read[j] = "N"
        self.rev_read = "".join(read)
        return(self.rev_read)

    def add_error_fwd(self):
        read = list(self.fwd_read)
        phred = list(self.fwd_phred)
        for j in range(0, len(read)):
            if npr.random() < quality.qdict[phred[j]]:
                read[j] = npr.choice(["A", "T", "G", "C"])
        self.fwd_err = "".join(read)
        return(self.fwd_err)

    def add_error_rev(self):
        read = list(self.rev_read)
        phred = list(self.rev_phred)
        for j in range(0, len(read)):
            if npr.random() < quality.qdict[phred[j]]:
                read[j] = npr.choice(["A", "T", "G", "C"])
        self.rev_err = "".join(read)
        return(self.rev_err)

    def combine_fwd(self):
        self.fwd_fq = f"@{self.name}_{self.id}/1\n{self.fwd_err}\n+\n{self.fwd_phred}\n"
        return(self.fwd_fq)

    def combine_rev(self):
        self.rev_fq = f"@{self.name}_{self.id}/2\n{self.rev_err}\n+\n{self.rev_phred}\n"
        return(self.rev_fq)


def prepare_run(all_frag, all_fwd_err, all_rev_err):
    res = [[i, j, k, l]
           for i, j, k, l in zip(all_frag, all_fwd_err, all_rev_err, range(0, len(all_frag)))]
    return(res)


def generate_fq(iterables, name,  mutate, mutrate, damage, geom_p, themin, themax, fwd_adaptor, rev_adaptor, read_length):
    frag = iterables[0]
    phred_fwd = iterables[1]
    phred_rev = iterables[2]
    num = iterables[3]
    frg = fragment(sequence=frag, num=num, name=name,
                   fwd_phred=phred_fwd, rev_phred=phred_rev)
    frg.reverse_complement()
    if mutate:
        frg.mutate_fwd(mutrate=mutrate)
        frg.mutate_rev(mutrate=mutrate)
    if damage:
        frg.add_damage_fwd(geom_p=geom_p, scale_min=themin, scale_max=themax)
        frg.add_damage_rev(geom_p=geom_p, scale_min=themin, scale_max=themax)
    frg.complement_read_fwd(adaptor=fwd_adaptor, read_length=read_length)
    frg.complement_read_rev(adaptor=rev_adaptor, read_length=read_length)
    frg.add_error_fwd()
    frg.add_error_rev()
    frg.combine_fwd()
    frg.combine_rev()
    return([frg.fwd_fq, frg.rev_fq])

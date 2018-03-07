#!/usr/bin/env python

from numpy import random as npr

QUALITY = "G"

def get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)

def reverse_complement(dna) :
    dna = dna.upper()
    '''
    Reverse complement a DNA string
    '''
    dna = dna[::-1]
    revcom = []
    complement = {"A" : "T", "T" : "A" , "G" : "C" , "C" : "G"}
    for letter in dna :
        for key in complement.keys() :
            if letter == key :
                revcom.append(complement[key])

    return "".join(revcom)

def read_fasta (file_name):
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
                result = result+line
    return([result, len(result)])

def random_insert(read_fasta_out, insert_lengths, read_length, minlen):
    genome = read_fasta_out[0]
    genome_length = read_fasta_out[1]
    result = []
    for i in insert_lengths:
        if i >= minlen:
            insert_start = npr.randint(0, genome_length-read_length)
            insert_end = insert_start + i + 1
            insert = genome[insert_start:insert_end]
            result.append(insert)
    return(result)

def complement_read(all_inserts, adaptor, read_length):
    result = []
    for insert in all_inserts:
        inlen = len(insert)
        if inlen < read_length:
            diff = read_length - inlen
            to_add = adaptor[0:diff]
            read = insert+to_add
        elif inlen == read_length:
            read = insert
        elif inlen > read_length:
            read = insert[0:read_length]
        result.append(read)
    return(result)

def add_error(all_reads, error_rate):
    for i in range(0, len(all_reads)):
        read = list(all_reads[i])
        for j in range(0, len(read)):
            if npr.random() < error_rate:
                read[j] = npr.choice(["A","T","G","C"])
                all_reads[i] = "".join(read)
    return(all_reads)

def prepare_fastq(fastq_dict, fwd_reads, rev_reads, basename, read_length):
    fastq_dict[basename] = [[] for i in range(2)]
    cnt = 1
    for read1, read2 in zip(fwd_reads, rev_reads):
        towrite_fwd = "@"+basename+"_"+str(cnt)+"/1"+"\n"+read1+"\n+\n"+QUALITY*read_length+"\n"
        fastq_dict[basename][0].append(towrite_fwd)
        towrite_rev = "@"+basename+"_"+str(cnt)+"/2"+"\n"+read2+"\n+\n"+QUALITY*read_length+"\n"
        fastq_dict[basename][1].append(towrite_rev)
        cnt += 1
    return(fastq_dict)

def write_fastq_multi(fastq_dict, outputfile):
    with open(outputfile+".1.fastq","w") as f1:
        with open(outputfile+".2.fastq","w") as f2:
            for akey in fastq_dict.keys():
                for reads1 in fastq_dict[akey][0]:
                    f1.write(reads1)
                for reads2 in fastq_dict[akey][1]:
                    f2.write(reads2)

def write_fastq(all_reads, basename, orientation, read_length, outfile):
    if not outfile:
        with open(basename+"."+str(orientation)+".fastq", "w") as fw:
            for read in all_reads:
                fw.write("@"+basename+"\n")
                fw.write(read+"\n")
                fw.write("+\n")
                fw.write(QUALITY*read_length+"\n")
    else:
        with open(outfile+"."+str(orientation)+".fastq", "w") as fw:
            for read in all_reads:
                fw.write("@"+basename+"\n")
                fw.write(read+"\n")
                fw.write("+\n")
                fw.write(QUALITY*read_length+"\n")



def run_read_simulation(INFILE, NREAD, COV, READLEN, INSERLEN, LENDEV, A1, A2, OUTFILE, MINLENGTH, ERR):
    print("INFILE: ", INFILE)
    if COV:
        print("COV: ", COV)
    else:
        print("NREAD: ", NREAD)
    print("READLEN: ", READLEN)
    print("INSERLEN: ", INSERLEN)
    print("LENDEV: ", LENDEV)
    print("A1: ", A1)
    print("A2: ", A2)
    print("OUTFILE: ", OUTFILE)

    nread = None


    basename = get_basename(INFILE)
    fasta = read_fasta(INFILE)

    if COV:
        nread = int(fasta[1]/INSERLEN)
        print("nread: ", nread)

    insert_lengths = [int(n) for n in npr.normal(INSERLEN, LENDEV, nread)]




    all_inserts = random_insert(fasta, insert_lengths, READLEN, MINLENGTH)
    fwd_inserts = all_inserts
    rev_inserts = [reverse_complement(i) for i in all_inserts]
    fwd_reads = complement_read(fwd_inserts, A1, READLEN)
    fwd_reads = add_error(fwd_reads, ERR)
    rev_reads = complement_read(rev_inserts, A2, READLEN)
    rev_reads = add_error(rev_reads, ERR)

    write_fastq(fwd_reads, basename, 1, READLEN, OUTFILE)
    write_fastq(rev_reads, basename, 2, READLEN, OUTFILE)

def run_read_simulation_multi(INFILE, NREAD, COV, READLEN, INSERLEN, LENDEV, A1, A2, MINLENGTH, ERR, fastq_dict):
    print("INFILE: ", INFILE)
    if COV:
        print("COV: ", COV)
    else:
        print("NREAD: ", NREAD)
    print("READLEN: ", READLEN)
    print("INSERLEN: ", INSERLEN)
    print("LENDEV: ", LENDEV)
    print("A1: ", A1)
    print("A2: ", A2)
    nread = None


    basename = get_basename(INFILE)
    fasta = read_fasta(INFILE)

    if COV:
        nread = int(fasta[1]/INSERLEN)
        print("nread: ", nread)

    insert_lengths = [int(n) for n in npr.normal(INSERLEN, LENDEV, nread)]




    all_inserts = random_insert(fasta, insert_lengths, READLEN, MINLENGTH)
    fwd_inserts = all_inserts
    rev_inserts = [reverse_complement(i) for i in all_inserts]
    fwd_reads = complement_read(fwd_inserts, A1, READLEN)
    fwd_reads = add_error(fwd_reads, ERR)
    rev_reads = complement_read(rev_inserts, A2, READLEN)
    rev_reads = add_error(rev_reads, ERR)

    prepare_fastq(fastq_dict = fastq_dict, fwd_reads = fwd_reads, rev_reads = rev_reads, basename = basename, read_length = READLEN)
    return(nread * INSERLEN)

def write_stat(stat_dict, stat_out):
    nbases = []
    for akey in stat_dict:
        nbases.append(stat_dict[akey])
    totbases = sum(nbases)
    with open(stat_out,"w") as fs:
        fs.write("Organism, percentage of metagenome\n")
        for akey in stat_dict:
            fs.write(akey+","+str(stat_dict[akey]/totbases)+"\n")

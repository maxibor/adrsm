import random
order = 10
readsize = 76


def generate_kmer(qualities, order, readsize):
    starts = []
    kmers = {}
    for q in qualities:
        for i in range(0, len(q)-order):
            kmer = q[i:i+order]
            if i == 0:
                starts.append(kmer)
            if kmer not in kmers.keys():
                kmers[kmer] = []
            kmers[kmer].append(q[i+order])
    return([starts, kmers])


def mchain(starts, kmers, readsize, order):
    seed = random.choice(starts)
    while(seed not in kmers.keys()):
        seed = random.choice(starts)
    chain = seed
    while(len(chain) < readsize):
        if seed in kmers.keys():
            first = kmers[seed]
            nxt = random.choice(first)
            chain += nxt
            seed = chain[len(chain)-order:len(chain)]

        else:
            return(0)
    return(chain)

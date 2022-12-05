#!/usr/bin/env python
import sys
import argparse
from concurrent import futures
from collections import namedtuple

import RNA
from Bio import SeqIO
from Bio.SeqUtils import GC

Probe = namedtuple('Probe', ['mfe', 'gc', 'seq', 'ss'])

def predictRNAFold(kmers):
    probes = list()
    for kmer in kmers:
        pred = RNA.fold(kmer.seq)
        probes.append(Probe(round(pred[1], 2), kmer.gc, kmer.seq, pred[0]))

    return probes

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def flatten(l):
    return [item for sublist in l for item in sublist]

def main():
    userInput = argparse.ArgumentParser("Find ISH Probes")
    userInput.add_argument('-f', '--file', action='store', required=True,
                            help='The FASTA file to find probes in')
    userInput.add_argument('-L', '--max-length', type=int, default=200,
                            help="Max length of probes")
    userInput.add_argument('-l', '--min-length', type=int, default=150,
                            help="Min length of probes.")
    userInput.add_argument('-p', '--process', type=int, default=1,
                            help="Number of process to folding RNA.")
    userInput.add_argument('-m', '--min-mfe', type=int, default=None,
                            help="The minimum allowed minimum free energy (MFE)")
    userInput.add_argument('-g', '--min_GC', dest = "min_gc", default=None, type=int,
                           help='The minimum allowed percent G + C')
    userInput.add_argument('-G', '--max_GC', dest = "max_gc", default=None, type=int,
                           help='The maximum allowed percent G + C')
    args = userInput.parse_args()

    for seq_record in SeqIO.parse(args.file, 'fasta'):
        # Get anti-sense strand
        block = str(seq_record.seq.reverse_complement()).upper()

        # Generage all possible k-mers.
        kmers = list()
        for k in range(args.min_length, args.max_length + 1):
            kmers += [block[i:i+k] for i in range(0, len(block) - k + 1)]

        #TODO: filter kmer with given GC limits.
        probes = [Probe(None, round(GC(kmer), 2), kmer, None) for kmer in kmers]

        if not args.min_gc is None:
            probes = list(filter(lambda x: x.gc > args.min_gc, probes))
        if not args.max_gc is None:
            probes = list(filter(lambda x: x.gc < args.max_gc, probes))

        # Parallel computing
        with futures.ProcessPoolExecutor(args.process) as executor:
            res = executor.map(predictRNAFold, list(chunks(probes, args.process)))

        # Combine and sort results
        probes = flatten(res)
        if not args.min_mfe is None:
            probes = list(filter(lambda x: x.mfe > args.min_mfe, probes))
        probes.sort(key=lambda x: x.mfe, reverse=True)

        # Output
        for probe in probes:
            sys.stdout.write("{seqid},{mfe},{gc},{kmer},{ss}\n".format(
                mfe = probe.mfe, gc = probe.gc, kmer = probe.seq, ss = probe.ss,
                seqid = seq_record.id
            ))


if __name__ == '__main__':
    main()

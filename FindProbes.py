#!/usr/bin/env python
import sys
import argparse
from concurrent import futures

import RNA
from Bio import SeqIO

def predictRNAFold(kmers):
    return [RNA.fold(kmer) + [kmer,] for kmer in kmers]

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
    userInput.add_argument('-g', '--min_GC', default=20, type=int,
                           help='The minimum allowed percent G + C, default is 20')
    userInput.add_argument('-G', '--max_GC', default=80, type=int,
                           help='The maximum allowed percent  G + C, default is 80')
    args = userInput.parse_args()

    for seq_record in SeqIO.parse(args.file, 'fasta'):
        block = str(seq_record.seq.reverse_complement()).upper()
        kmers = list()
        for k in range(args.min_length, args.max_length):
            kmers += [block[i:i+k] for i in range(0, len(block) - k + 1)]
        with futures.ProcessPoolExecutor(args.process) as executor:
            res = executor.map(predictRNAFold, list(chunks(kmers, args.process)))
        predicts = flatten(res)
        predicts.sort(key=lambda x: x[1], reverse=True)

        for pred in predicts:
            sys.stdout.write("{seqid},{mfe},{kmer},{ss}\n".format(
                mfe = pred[1], kmer = pred[2], ss = pred[0],
                seqid = seq_record.id
            ))


if __name__ == '__main__':
    main()

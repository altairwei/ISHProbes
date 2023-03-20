#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
import csv
from concurrent import futures
from collections import namedtuple

import RNA
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp

Probe = namedtuple('Probe', ['mfe', 'gc', 'tm', 'seq', 'AB', 'AA', 'BB', 'A', 'B'])
NNTABLE = MeltingTemp.RNA_NN3

def runRNACofold(probes):
    # Run RNAcofold for antisense and sense probes
    cofold_inputs = [
        str(probe.seq) + "&" + str(Seq(probe.seq).reverse_complement()) for probe in probes]
    result = subprocess.check_output(
        ["RNAcofold", "--jobs=1",
         "-a", "-d2", "--noLP", "--noPS",
         "--output-format=D"],
        input="\n".join(cofold_inputs).encode())
    reader = csv.reader(result.decode().splitlines())
    next(reader, None) # Skip header
    cofolds = list(reader)
    # Extract results from RNAcofold
    outputs = list()
    for probe, cofold in zip(probes, cofolds):
        outputs.append(Probe(
            probe.mfe, probe.gc, probe.tm, probe.seq,
            # Note: Two columns were absent in the CSV header of RNAcofold v2.5.1,
            # so we only use integer to index this table.
            round(float(cofold[9]), 2), round(float(cofold[10]), 2),
            round(float(cofold[11]), 2), round(float(cofold[12]), 2),
            round(float(cofold[13]), 2)
        ))
    
    return outputs

def predictRNAFold(kmers):
    probes = list()
    for kmer in kmers:
        pred = RNA.fold(kmer.seq)
        probes.append(Probe(round(pred[1], 2), kmer.gc, kmer.tm, kmer.seq,
                            0, 0, 0, 0, 0))
    probes = runRNACofold(probes)
    return probes

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def flatten(lst):
    return [item for sublist in lst for item in sublist]

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
    userInput.add_argument('-H', '--header', default=False, action="store_true",
                            help='Write a header line to csv outputs.')
    args = userInput.parse_args()

    if args.header:
        sys.stdout.write("Seq_ID,MFE,GC,TM,AB_Heterodimer,AA_Homodimer,BB_Homodimer,A_Monomer,B_Monomer,K-mer\n")

    for seq_record in SeqIO.parse(args.file, 'fasta'):
        block = str(seq_record.seq).upper()

        # Generage all possible k-mers.
        kmers = list()
        for k in range(args.min_length, args.max_length + 1):
            kmers += [block[i:i+k] for i in range(0, len(block) - k + 1)]

        # Nearest neighbor thermodynamic Tm predictions for DNA do not really
        # extend beyond duplexes of ~60 bp, so the Tm values at 120 would not
        # be accurate.
        probes = [Probe(
            None, round(gc_fraction(kmer) * 100, 2),
            round(MeltingTemp.Tm_NN(kmer, nn_table = NNTABLE), 2),
            kmer, 0, 0, 0, 0, 0) for kmer in kmers]

        if not args.min_gc is None:
            probes = list(filter(lambda x: x.gc > args.min_gc, probes))
        if not args.max_gc is None:
            probes = list(filter(lambda x: x.gc < args.max_gc, probes))

        # Parallel computing
        if args.process > 1:
            with futures.ProcessPoolExecutor(args.process) as executor:
                res = executor.map(predictRNAFold, list(chunks(probes, args.process)))
        else:
            res = predictRNAFold(probes)

        # Combine and sort results
        probes = flatten(res)
        if not args.min_mfe is None:
            probes = list(filter(lambda x: x.mfe > args.min_mfe, probes))
        probes.sort(key=lambda x: x.mfe, reverse=True)

        # Output
        for probe in probes:
            sys.stdout.write("{seqid},{mfe},{gc},{tm},{AB},{AA},{BB},{A},{B},{kmer}\n".format(
                mfe = probe.mfe, gc = probe.gc, kmer = probe.seq,
                seqid = seq_record.id, tm = probe.tm,
                AB = probe.AB, AA = probe.AA, BB = probe.BB, A = probe.A, B = probe.B
            ))

    # Clean *dot5.ps files.
    for psfile in ["AAdot5.ps", "ABdot5.ps", "Adot5.ps",
                   "BBdot5.ps", "Bdot5.ps"]:
        if os.path.isfile (psfile):
            os.remove (psfile)


if __name__ == '__main__':
    main()

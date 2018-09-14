import Fasta_one_line as fol
import multiprocessing as mp
import os
import argparse
from subprocess import Popen


def make_fasta(seqName, seqDict):
    locus = seqName.split("|")[1].split(".")[0]
    fastaName = "temp_plasmid_seq.{}.fna".format(locus)
    with open(fastaName, 'w') as f:
        f.write(">{}\n{}".format(seqName, seqDict[seqName]))
    return fastaName


def einverted_cmd(fasta, gapLength):
    num = fasta.split(".")[1]
    outname = "temp_plasmid_inverted.{}".format(num)
    outfasta = outname + ".fasta"
    cmd = "/work/software/EMBOSS-6.6.0/emboss/einverted -sequence {} -outfile {} -gap 12 -threshold 50 -match 3 " \
          "-mismatch -4 -outseq {} -maxrepeat {}".format(fasta, outname, outfasta, gapLength)
    proc = Popen(cmd, shell=True)
    proc.wait()
    return outname


parser = argparse.ArgumentParser(description="Apply einverted in parallel to allow large gap sizes")

parser.add_argument("Fasta", help="Name of fasta file to run einverted on")
parser.add_argument("-g", "--gap", default=10000, type=int, help="Desired max gap length")
parser.add_argument("-p", "--processes", default=2, type=int, help="Number of processes to use")

args = parser.parse_args()

fastaDict = fol.one_line_d(args.Fasta)
fastaList = list(fastaDict.keys())

with mp.Pool(processes=args.processes) as pool:
    fastas = [pool.apply_async(make_fasta, args=(name, fastaDict)) for name in fastaList]
    fastaFileList = [p.get() for p in fastas]
    einverteds = [pool.apply_async(einverted_cmd, args=(name, args.gap)) for name in fastaFileList]
    invertedOut = [p.get() for p in einverteds]

map(os.remove, fastaFileList)

cat_cmd = "cat temp_plasmid_inverted.* > plasmid_inverted_parallel"
cat = Popen(cat_cmd, shell=True)
cat.wait()
map(os.remove, invertedOut)

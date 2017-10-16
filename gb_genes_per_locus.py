from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Count coding regions per locus in GenBank file")

parser.add_argument("GenBankFile", help="Input file name")
parser.add_argument("Outfile", help="Output file name")

args = parser.parse_args()

records = SeqIO.parse(args.GenBankFile, 'genbank')
outfile = open(args.Outfile, 'w')

outfile.write("locus\tCDS\n")
for record in records:
    plasmidName = record.name

    codingRegions = 0

    for cds in record.features:
        if cds.type == "CDS":
            codingRegions += 1

    outfile.write("{}\t{}\n".format(plasmidName, codingRegions))
outfile.close()
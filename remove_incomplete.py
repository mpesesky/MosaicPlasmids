from Bio import SeqIO
import Fasta_one_line as Fol
import argparse

parser = argparse.ArgumentParser(description="Remove plasmids by accession prefix")

parser.add_argument("Infile")
parser.add_argument("Outfile")

args = parser.parse_args()


if args.Infile.endswith("fna"):
    outfile = open(args.Outfile, 'w')
    seqs = Fol.one_line_d(args.Infile)
    for header in seqs.keys():
        if ("NC_" in header) or ("NZ_" in header):
            outfile.write(">{}\n{}\n".format(header, seqs[header]))
    outfile.close()
elif args.Infile.endswith("gbff"):
    outlist = []
    records = SeqIO.parse(args.Infile, 'genbank')
    for record in records:
        if record.id.startswith("NC_") or record.id.startswith("NZ_"):
            outlist.append(record)
    SeqIO.write(outlist, args.Outfile, "genbank")
else:
    print("Extension not recognized")
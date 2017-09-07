import argparse
import Fasta_one_line as fol

def query_to_locus(query):
    return query.split("|")[1]

def header_to_query(header):
    return header.split(" ")[0].lstrip(">")

parser = argparse.ArgumentParser(description="Split a set of sequences based on their presence in a BLAST file")

parser.add_argument("FastaFile", help="File of sequences to be split")
parser.add_argument("BlastFile", help="File with subset of Fasta")
parser.add_argument("-m", "--matches", type=str, default=None, help="Matching output")
parser.add_argument("-n", "--non_matches", type=str, default=None, help="Non-Match output")
parser.add_argument("-c", "--combo", type=str, default="Output.txt", help="Combined output (default = Output.txt)")
args = parser.parse_args()

headers = list(map(header_to_query, fol.one_line_d(args.FastaFile).keys()))
#print(headers[:5])

blastFile = open(args.BlastFile)
matches = []

for line in blastFile:
    if line.startswith("qseqid"):
        continue
    fields = line.split("\t")
    query = fields[0]
    subject = fields[1]
#    print(query)
#    exit()

    if subject in headers:
        matches.append(query_to_locus(subject))
        headers.remove(subject)
blastFile.close()

nonMatches = map(query_to_locus, headers)

if args.matches is not None:
    outMatch = open(args.matches, 'w')
    for item in matches:
        outMatch.write(item)
        outMatch.write("\n")
    outMatch.close()

    outNon = open(args.non_matches, 'w')
    for item in nonMatches:
        outNon.write(item)
        outNon.write("\n")
    outNon.close()
else:
    output = open(args.combo, 'w')
    for item in matches:
        output.write("{}\t{}\n".format(item, "Present"))
    for item in nonMatches:
        output.write("{}\t{}\n".format(item, "Absent"))
    output.close()

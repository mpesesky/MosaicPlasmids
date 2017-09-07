import Fasta_one_line as fol

seqDict = fol.one_line_d("ncbi_plasmid.fna")
outfile = open("Plasmid_lengths.txt", 'w')
outfile.write("locus\tlength\n")

for header in seqDict.keys():
    locus = header.split("|")[1].split(".")[0]
    outfile.write("{}\t{}\n".format(locus,len(seqDict[header])))
outfile.close()


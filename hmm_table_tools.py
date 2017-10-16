import pandas as pd
import Fasta_one_line as fol
from Bio import SeqIO

colNames = ['target name', 'taccession', 'query name', 'qaccession', 'fE-value', 'fscore', 'fbias', 'dE-value',
            'dscore', 'dbias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description of target']

def import_hmm_table(tableFile):
    df = pd.read_fwf(tableFile, names=colNames, comment='#')
    return df


def min_eval(hmmDF):
    outdf = pd.DataFrame(columns=colNames)

    for index, match in hmmDF.iterrows():
        if match['query name'] not in outdf['query name'].values:
            outdf = outdf.append(match)
    return outdf


def filter_eval(hmmDF, maxEval):
    outdf = hmmDF[hmmDF['fE-value'].map(float) <= maxEval]
    return outdf


def import_fasta(fastaFile):
    seqDict = fol.blast_dict(fastaFile)
    longList = fol.one_line_d(fastaFile).keys()
    descDict = {}
    lenDict = {}
    for header in seqDict.keys():
        lenDict[header] = len(seqDict[header])
    for header in longList:
        parts = header.rstrip("\n").lstrip(">").split(" ")
        descDict[parts[0]] = " ".join(parts[1:])
    return lenDict, descDict, seqDict


def add_lengths(hmmDF, lengths):
    hmmDF['qlen'] = hmmDF['query name'].map(lengths)
    return hmmDF


def add_descriptions(hmmDF, desc):
    hmmDF['description of target'] = hmmDF['query name'].map(desc)
    return hmmDF


def make_tsv(tableFile, tsvName):
    if tsvName[-4:] != ".tsv":
        tsvName = tsvName + ".tsv"
    print("reprocessing the hmmscan table to a tab separated value file that is clean")
    with open(tableFile, 'rt') as input, open(tsvName, 'wt') as output:
        for line in input:
            fields = line.split()
            if fields:
                if fields[0][0] != '#':
                    if len(fields) < 19:
                        raise Exception("Incorrect number of fields",
                                        line)
                    else:
                        for field in fields[:18]:
                            output.write("%s\t" % field.replace('"', ''))
                        for field in fields[18:len(fields) - 1]:
                            output.write("%s " % field.replace('"', ''))
                        output.write("%s\n" % fields[len(fields) - 1].replace('"', ''))
                else:
                    pass
    return tsvName


def import_tsv(tableFile):
    df = pd.read_table(tableFile, sep="\t", names=colNames, comment='#')
    return df


def gb_add(hmmDF, gbFile):
    plasmids = SeqIO.parse(gbFile, 'genbank')
    proteins = hmmDF['query name'].tolist()
    hmmDF['plasmid ID'] = "i"
    outDF = pd.DataFrame(data=None, columns=hmmDF.columns)

    for plasmid in plasmids:
        plasid = plasmid.name
        for gene in plasmid.features:
            if gene.type == 'CDS':
                try:
                    gene_id = "ref|{}|".format(gene.qualifiers['protein_id'][0])
                except KeyError:
                    continue
                if gene_id in proteins:
                    geneRow = hmmDF[hmmDF['query name'] == gene_id]
                    geneRow['plasmid ID'] = plasid
                    outDF = outDF.append(geneRow)
    return outDF


def extract_seqs(hmmDF, seqDict):
    retDict ={}

    targetList = list(hmmDF['target name'])

    for header in seqDict.keys():
        if header in targetList:
            retDict[header] = seqDict[header]
    return retDict


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Manipulate HMMER table output")

    parser.add_argument("TableFile", help="HMMER Table output file")
    parser.add_argument("-e", "--min_eval", metavar="OUT", default=None, type=str, help="Print table to OUT with only "
                                                                                        "minimum e-values per query to OUT")
    parser.add_argument("-f", "--fasta", default=None, type=str, help="Import Fasta for length information")
    parser.add_argument("-m", "--filter_eval", default=1, type=float, help="Set minimum e-value")
    parser.add_argument("-g", "--gbfile", type=str, default=None, help="GenBank file to match proteins to plasmids")
    parser.add_argument("-p", "--protein", type=str, default=None, help="Export Matching seqs to PROTEIN file")

    args = parser.parse_args()

    if os.path.splitext(args.TableFile)[1] == ".txt":
        baseName = os.path.splitext(args.TableFile)[0]
        args.TableFile = make_tsv(args.TableFile, baseName)

    hmm = import_tsv(args.TableFile)

    if args.filter_eval < 1.0:
        hmm = filter_eval(hmm, args.filter_eval)

    if args.fasta is not None:
        lengths, descriptions, sequences = import_fasta(args.fasta)
        hmm = add_lengths(hmm, lengths)
        hmm = add_descriptions(hmm, descriptions)

    if args.min_eval is not None:
        outdf = min_eval(hmm)
        print(outdf)
        if args.gbfile is not None:
            outdf = gb_add(outdf, args.gbfile)
        outdf.to_csv(args.min_eval, sep="\t")

    if args.protein is not None:
        try:
            outSeqs = extract_seqs(hmm, sequences)
        except ValueError:
            print("Must give Fasta to get Fasta")
            exit()
        outfile = open(args.protein, 'w')
        for header in sequences.keys():
            outfile.write(">{}\n{}\n".format(header, sequences[header]))

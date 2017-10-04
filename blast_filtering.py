import pandas as pd
import math
import Fasta_one_line as fol

colNames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                "send", "evalue", "bitscore"]

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', 15)
    print(x.to_csv(sep="\t", index=False))
    pd.reset_option('display.max_rows')
    pd.reset_option('expand_frame_repr')

def lengthFilter(valueTable, minLen):
    newdf = valueTable[valueTable['length'].map(float) > minLen]
    return newdf

def nonself(valueTable):
    newdf = valueTable[valueTable['qseqid'] != valueTable['sseqid']]
    return newdf

def info(valueTable):
    newdf = valueTable.mode()
    newdf.append(valueTable.mean(numeric_only=True), ignore_index=True)
    newdf.append(valueTable.max(numeric_only=True), ignore_index=True)
    return newdf

def abundance(valueTable, numbins):
    queryCounts = valueTable['qseqid'].value_counts()
    binsize = int(math.ceil(float(queryCounts.max())/numbins))

    print(queryCounts[queryCounts > 1200].index.values)

    abundanceDict = {}
    lastAmount = 0
    for i in range(1,(numbins+1)):
        currentBin = binsize * i
        currentAmount = queryCounts[queryCounts <= currentBin].count()
        abundanceDict[currentBin] = currentAmount - lastAmount
        lastAmount = currentAmount
    x = pd.Series(abundanceDict).sort_index()
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')
    return x

def queryLengths(seqFile):
    seqDict = fol.one_line_d(seqFile)
    lenDict = {}
    for header in seqDict.keys():
        query = header.split(" ")[0].lstrip(">")
        lenDict[query] = len(seqDict[header])
    return lenDict

def percFilter(valueTable, lenDict, perc):
    valueTable['qlen'] = valueTable['qseqid'].map(lenDict)
    valueTable['slen'] = valueTable['sseqid'].map(lenDict)
    retTable = valueTable[((valueTable['length'].map(float) / valueTable['qlen'].map(float)) >= (perc / 100))|((valueTable['length'].map(float) / valueTable['slen'].map(float)) >= (perc / 100))]
    return retTable[colNames]

def percidFilter(valueTable, percid):
    newdf = valueTable[valueTable['pident'].map(float) > percid]
    return newdf

def otherBLASTFilter(mainTable, filterTable):
    keys = ['qseqid', 'sseqid']
    i1 = mainTable.set_index(keys).index
    i2 = filterTable.set_index(keys).index
    return mainTable[~i1.isin(i2)]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter default tabular BLAST results. Prints to command line.")

    parser.add_argument("-p", "--prot", action='store_true', help="Results are blastp not blastn")
    parser.add_argument("-n", "--nonself", action='store_true', help="Remove self matches")
    parser.add_argument("-l", "--length", default=None, help="Remove matches under LENGTH bp long")
    parser.add_argument("-i", "--info", action='store_true', help="Only output useful information on BLAST result file")
    parser.add_argument("BLAST", help="BLAST results generated with outfmt=6")
    parser.add_argument("-q", "--queryAb", metavar="X", default=0, type=int, help="Chart query abundance into X even bins")
    parser.add_argument("-c", "--percent", type=float, help="Filter out matches below percent length")
    parser.add_argument("-s", "--seq", type=str, default=None, help="Query seq file, to determine query lengths")
    parser.add_argument("-e", "--percid", type=float, default=None, help="Filter out matches lower than percid")
    parser.add_argument("-o", "--other", type=str, default=None, help="Filter out query subject pairs present in second BLAST file")

    args = parser.parse_args()

    results = pd.read_table(args.BLAST, sep="\t", header=None, names=colNames)

    if args.queryAb > 0:
        outframe = abundance(results, args.queryAb)
        exit()

    if args.info:
        information = info(results)
        print_full(information)
        print(results['qseqid'].unique().size)
        exit()

    if args.nonself:
        results = nonself(results)

    if args.percid is not None:
        results = percidFilter(results, args.percid)

    if args.length is not None:
        results = lengthFilter(results, int(args.length))

    if args.percent is not None:
        if args.seq is None:
            print("Need seq file to determine percent query match")
        else:
            lengths = queryLengths(args.seq)
            results = percFilter(results, lengths, args.percent)

    if args.other is not None:
        otherBlast = pd.read_table(args.other, sep="\t", header=None, names=colNames)
        results = otherBLASTFilter(results, otherBlast)

    print_full(results)
import pandas as pd
from Bio import SeqIO
import argparse
import re

class Plasmid:
    def __init__(self, acc):
        self.accession = acc
        self.first_starts = []
        self.first_ends = []
        self.second_starts = []
        self.second_ends = []

    def __str__(self):
        return "{}\nFirst starts: {}\nFirst Ends: {}\nSecond Starts: {}\nSecond Ends: " \
               "{}\n".format(self.accession, self.first_starts, self.first_ends, self.second_starts, self.second_ends)

def import_repeats(aln):
    plasDict = {}

    seqLine = re.compile(r"(?P<start>\d+)\s\w+\s(?P<end>\d+)")
    defLine = re.compile(r"(?P<accession>.+):\sScore")

    with open(aln) as o:
        curplas = ""
        first = False
        for line in o:
            isDefLine = defLine.search(line)

            if isDefLine is not None:
                curplas = isDefLine.groupdict()['accession']
                first = True
                if curplas not in plasDict:
                    plasDict[curplas] = Plasmid(curplas)
            else:
                isSeqLine = seqLine.search(line)

                if isSeqLine is not None:
                    nums = isSeqLine.groupdict()
                    if first:
                        first = False
                        plasDict[curplas].first_starts.append(int(nums['start']))
                        plasDict[curplas].first_ends.append(int(nums['end']))
                    else:
                        plasDict[curplas].second_starts.append(int(nums['start']))
                        plasDict[curplas].second_ends.append(int(nums['end']))

    return plasDict


def seqid_to_locus(seqid):
    return seqid.split("|")[1]


def seqid_to_short_locus(seqid):
    return seqid.split("|")[1].split(".")[0]


def match_repeats(row, repeatDF, dist, start):
    try:
        plasmid = repeatDF[row['qlocus']]
    except KeyError:
        return 'NA'

    if start:
        pos = int(row['qstart'])
        repList1 = plasmid.first_starts
        repList2 = plasmid.first_ends
    else:
        pos = int(row['qend'])
        repList1 = plasmid.second_ends
        repList2 = plasmid.second_starts

    for i in range(len(repList1)):
        repStart = repList1[i]
        repEnd = repList2[i]
        if ((pos >= repStart) and (pos <= repEnd)) or (abs(pos - repStart) <= dist) or (abs(pos - repEnd) <= dist):
            return i
    return 'NA'


def add_pos(row, gbDict, start):
    if row['plasmid ID'].startswith("NT") or row['plasmid ID'].startswith("NG"):
        return 'NA'
    try:
        record = gbDict["{}.1".format(row['plasmid ID'])]
    except KeyError:
        try:
            record = gbDict["{}.2".format(row['plasmid ID'])]
        except KeyError:
            try:
                record = gbDict["{}.3".format(row['plasmid ID'])]
            except KeyError:
                try:
                    record = gbDict["{}.4".format(row['plasmid ID'])]
                except KeyError:
                    print(row['plasmid ID'])
                    exit()
    for feat in record.features:
        if feat.type != 'CDS':
            continue
        if ('pseudo' in feat.qualifiers) or ('pseudogene' in feat.qualifiers):
            continue
        try:
            featid = feat.qualifiers['protein_id'][0]
        except KeyError:
            print(feat.qualifiers)
            exit()
        if featid == row['qlocus']:
            if start:
                return feat.location.start
            else:
                return feat.location.end
    return 'NA'


def between(x, first, last):
    return (x >= first) and (x <= last)


def find_transposases(row, transDF):
    plasmidName = row['qlocus'].split(".")[0]
    releDF = transDF[transDF['plasmid ID'] == plasmidName]
    starts = releDF['qstart'].tolist()
    ends = releDF['qend'].tolist()
    count = 0
    for i in range(len(starts)):
        if between(starts[i], row['qstart'], row['qend']) or between(ends[i], row['qstart'], row['qend']):
            count += 1
    return count


parser = argparse.ArgumentParser(description="Find inverted repeats near mosaic edges")

parser.add_argument("Filtered", help="BLAST 6 output delineating mosaic fragments")
parser.add_argument("Inverted", help="Alignment of einverted")
parser.add_argument("Outfile", help="Desired output table name")
parser.add_argument("-d", "--dist", default=5, type=int, help="Distance from repeat to find edges")
parser.add_argument("-t", "--transposases", default=None, help="HMMER output of transposases locations")
parser.add_argument("-g", "--genbank", default=None, help="Genbank file of plasmids")

args = parser.parse_args()

repeats = import_repeats(args.Inverted)
#test = list(repeats.keys())[0]
#print(repeats[test])

fragments = pd.read_table(args.Filtered, sep="\t")

fragments['qlocus'] = fragments['qseqid'].map(seqid_to_locus)
fragments['slocus'] = fragments['sseqid'].map(seqid_to_locus)

fragments['Start_repeat'] = fragments.apply(lambda x: match_repeats(x, repeats, args.dist, True), axis=1)

fragments['End_repeat'] = fragments.apply(lambda x: match_repeats(x, repeats, args.dist, False), axis=1)

if (args.transposases is not None) and (args.genbank is not None):
    trans = pd.read_table(args.transposases, index_col=0, sep="\t")
    gb = SeqIO.to_dict(SeqIO.parse(args.genbank, 'genbank'))
    trans['qlocus'] = trans['query name'].map(seqid_to_locus)
    trans['qstart'] = trans.apply(lambda x: add_pos(x, gb, True), axis=1)
    trans['qend'] = trans.apply(lambda x: add_pos(x, gb, False), axis=1)

    fragments['Transposases'] = fragments.apply(lambda x: find_transposases(x, trans), axis=1)

    fragments[['qlocus', 'length', 'qstart', 'qend', 'slocus', 'sstart', 'send', 'Start_repeat', 'End_repeat',
               'Transposases']].to_csv(args.Outfile, sep="\t")
else:
    fragments[['qlocus', 'length', 'qstart', 'qend', 'slocus', 'sstart', 'send', 'Start_repeat',
               'End_repeat']].to_csv(args.Outfile, sep="\t")
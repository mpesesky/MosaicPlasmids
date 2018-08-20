import pandas as pd
import argparse


def seqid_to_locus(seqid):
    return seqid.split("|")[1].split(".")[0]

def seqid_to_rep(seqid):
    return seqid.split("_")[0].split("(")[0]

def row_overlap(row1, row2, cutoff):
    x = float(row1.send - row1.sstart)
    y = float(row2.send - row2.sstart)
    if y > x:
        x = y

    if (row1.sstart <= row2.sstart) and (row1.send >= row2.send):
        return True
    elif (row2.sstart <= row1.sstart) and (row2.send >= row1.send):
        return True
    elif (row1.sstart <= row2.sstart) and (row1.send > row2.sstart):
        if ((row1.send - row2.sstart)/x) >= cutoff:
            return True
        else:
            return False
    elif (row2.sstart <= row1.sstart) and (row2.send > row1.sstart):
        if ((row2.send - row1.sstart) / x) >= cutoff:
            return True
        else:
            return False
    else:
        return False


def remove_overlaps(grouped, cutoff):
    newDF = pd.DataFrame(columns=grouped.columns)

    for i, row in grouped.iterrows():
        if newDF.empty:
            newDF = newDF.append(row)
            continue

        duplicate = False
        for i2, row2 in newDF.iterrows():
            if row_overlap(row, row2, cutoff):
                duplicate = True
                if row['bitscore'] > row2['bitscore']:
                    newDF = newDF.drop(row2.name)
                    newDF = newDF.append(row)
                break

        if not duplicate:
            newDF = newDF.append(row)
    return newDF


parser = argparse.ArgumentParser(description="Create a column of rep genes for the MasterTable")

parser.add_argument("RepBlast", help="Filtered results of plasmid finder positive BLAST")
parser.add_argument("OutFile", help="Desired output file name")
parser.add_argument("-c", "--cutoff", default=0.8, help="Overlap cutoff value")

args = parser.parse_args()

blastDF = pd.read_table(args.RepBlast, sep="\t", dtype={'sstart': int, 'send': int, 'bitscore': float})

seqidList = list(blastDF['sseqid'].unique())

condensedDF = pd.DataFrame(columns=blastDF.columns)
for seqid in seqidList:
    condensedDF = condensedDF.append(remove_overlaps(blastDF[blastDF['sseqid'] == seqid], args.cutoff))

maxOnly = condensedDF.sort_values('bitscore').drop_duplicates(['sseqid'])[['sseqid', 'qseqid']]

repDF = condensedDF.groupby('sseqid')['qseqid'].count()

merged = pd.merge(pd.DataFrame(repDF), maxOnly, left_index=True, right_on='sseqid')
merged['locus'] = merged['sseqid'].map(seqid_to_locus)
merged['Best_rep_name'] = merged['qseqid_y'].map(seqid_to_rep)
merged.rename(columns={'qseqid_x': 'Rep_num'}, inplace=True)
merged[['locus', 'Rep_num', 'Best_rep_name']].to_csv(args.OutFile, sep="\t", index=False)

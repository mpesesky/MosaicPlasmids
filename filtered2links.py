#!/usr/bin/env python3

FILTER_LENGTH = 500

import csv

import pandas as pd

print("reading in contig names and lengths")
# read the file that contains lengths and names
#chr - NC_003894 NC_003894 1 9253 145,255,145
kdf = pd.read_csv("karyotype.txt", sep=' ', header=None, 
        names=["chr", "ignore", "contig", "name",
                "start", "end", "color"])
kdf["length"] = (kdf.end - kdf.start) + 1
kdf.drop(["chr", "ignore", "contig", "color", "start", "end"], 1, inplace=True)

print("read filtered matches")
# read the matches and fix up the names
mdf = pd.read_csv("../taxonAnalysis/mosaic_fragments_BLAST.txt",
        sep='\t')
fixed = mdf["qseqid"].apply(lambda x : x.split("|")[1].split(".")[0])
mdf["qseqid"] = fixed
fixed = mdf["sseqid"].apply(lambda x : x.split("|")[1].split(".")[0])
mdf["sseqid"] = fixed
mdf.rename(columns = { "length" : "alength" }, inplace=True)

print("reading in organism names")
odf = pd.read_csv("organisms.txt", sep='\t', header=None, 
        names=["name", "organism"])

print("merging data frames")
# join mdf twice, once with lengths on qseqid and then on sseqid
df = pd.merge(mdf, kdf, left_on="qseqid", right_on="name")
df.drop("name", 1, inplace=True);
df.rename(columns = { "length" : "qlength"}, inplace=True)
df = pd.merge(df, kdf, left_on="sseqid", right_on="name")
df.drop("name", 1, inplace=True);
df.rename(columns = { "length" : "slength"}, inplace=True)
# now join to organisms
df = pd.merge(df, odf, left_on="qseqid", right_on="name")
df.drop("name", 1, inplace=True);
df.rename(columns = { "organism" : "qorganism"}, inplace=True)
df = pd.merge(df, odf, left_on="sseqid", right_on="name")
df.drop("name", 1, inplace=True);
df.rename(columns = { "organism" : "sorganism"}, inplace=True)
df["color"] = df.apply(lambda x : "color=black" if x['qorganism'] == x['sorganism'] else "color=red", axis=1)

print(df.head())

# deprecated; all the processing is done by mitch's code
"""
print("filtering out redundant plasmids")
# remove any matches that are > 50% length of either
qssame = df["qseqid"] == df["sseqid"]
qscrit = df["alength"] < df["slength"] * .5
qqcrit = df["alength"] < df["qlength"] * .5
scrit = df["alength"] > 0
scrit = df["alength"] > 1000
scrit = df["alength"] > 5000
print(df.size)
df = df[(qssame) & (qscrit) & (qqcrit) & (scrit)]
print(df.size)
"""
# optional size handling is here
scrit = df["alength"] > FILTER_LENGTH
print(df.size)
df = df[(scrit)]
print(df.size)

print("writing to 'links.txt'")
links = df[["qseqid", "qstart", "qend", "sseqid", "sstart", "send", "color"]].copy()
links.to_csv("links.txt", sep='\t', quoting=csv.QUOTE_NONE, header=False, index=False)

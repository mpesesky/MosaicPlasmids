import pandas as pd
import argparse
import numpy as np
from multiprocessing import Pool


def inter_inc(row, master):
    try:
        query = master[master['locus'] == row['locus1']]['Best_rep_name'].tolist()[0]
        subject = master[master['locus'] == row['locus2']]['Best_rep_name'].tolist()[0]
    except IndexError:
        return 'color=blue'
    if pd.isnull(query) and pd.isnull(subject):
        return 'color=blue'

    if query == subject:
        return 'color=black'
    else:
        return 'color=red'


def run_inter_inc(df, otherDF):
    df['InterIncType'] = df.apply(lambda x: inter_inc(x, otherDF), axis=1)
    return df


def parallelize(df, otherDF, proc):
    data_split = np.array_split(df, proc)
    pool = Pool(proc)
    df = pd.concat(pool.map(lambda x: run_inter_inc(x, otherDF), data_split))
    pool.close()
    pool.join()
    return df


parser = argparse.ArgumentParser(description="Enumerate connections between inc groups")

parser.add_argument("Master", help="Master Plasmid Table")
parser.add_argument("Links", help="Mosaic link network table")
parser.add_argument("Output", help="Desired output file name")
parser.add_argument("-p", "--processes", default=1, type=int, help="Number of parallel processes to run")

args = parser.parse_args()

columns = ['locus1', 'start1', 'end1', 'locus2', 'start2', 'end2', 'Intergenus']

mdf = pd.read_table(args.Master, index_col=0, sep="\t")
ldf = pd.read_table(args.Links, sep="\t", names=columns)

#ldf = parallelize(ldf, mdf, args.processes)
ldf['InterIncType'] = ldf.apply(lambda x: inter_inc(x, mdf), axis=1)

ldf.drop(columns='Intergenus', inplace=True)

ldf.to_csv(args.Output, sep="\t", index=False, header=False)

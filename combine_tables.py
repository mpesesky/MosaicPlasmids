import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="Combine tab delimited text tables into one file")

parser.add_argument("Table1", help="First table file")
parser.add_argument("-t", "--tables", nargs="+", help="Additional Table Files")
parser.add_argument("-o", "--output", default="CombinedTable.txt", help="Output filename (default = CombinedTable.txt)")
parser.add_argument("-c", "--column", default="locus", help="Name of column to merge on (default = 'locus')")

args = parser.parse_args()

df1 = pd.read_csv(args.Table1, delimiter="\t")
if args.tables is not None:
    otherDFs = []
    for tableFile in args.tables:
        df1 = pd.merge(df1, pd.read_csv(tableFile, delimiter="\t"), how="outer", on=args.column)

df1.fillna('na', inplace=True)
df1.to_csv(args.output, sep="\t")
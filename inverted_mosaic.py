import pandas as pd
import argparse


def trans_count(dataframe):
    return len(dataframe[dataframe['Transposases'] > 0].index)


parser = argparse.ArgumentParser(description="Get stats on inverted repeats and mosaic fragments")

parser.add_argument("RepeatFile", help="File containing repeat and mosaic information")

args = parser.parse_args()

df = pd.read_table(args.RepeatFile, sep="\t", index_col=0)
denom = float(len(df.index))

hasRepeats = df[~(df['Start_repeat'].isnull()) & ~(df['End_repeat'].isnull())]
oneRepeat = df[(~(df['Start_repeat'].isnull()) & (df['End_repeat'].isnull())) |
               (df['Start_repeat'].isnull() & ~(df['End_repeat'].isnull()))]
noRepeats = df[(df['Start_repeat'].isnull()) & (df['End_repeat'].isnull())]
bracketed = hasRepeats[hasRepeats['Start_repeat'] == hasRepeats['End_repeat']]

print("Overall percent of mosaic fragments with at least one transposase: {}".format((trans_count(df)/denom)*100))

print("Percent of mosaic fragments with matching inverse repeats: {}".format((len(bracketed.index)/denom)*100))
print("And at least one transposase: {}".format((trans_count(bracketed)/denom)*100))

print("Percent of mosaic fragments with unmatched inverse repeats: "
      "{}".format(((len(hasRepeats.index) - len(bracketed.index))/denom)*100))
print("And at least one transposase: {}".format(((trans_count(hasRepeats) - trans_count(bracketed))/denom)*100))

print("Percent of mosaic fragments with only one inverted repeat: {}".format((len(oneRepeat.index)/denom)*100))
print("And at least one transposase: {}".format(((trans_count(oneRepeat) - trans_count(bracketed))/denom)*100))

print("Percent of mosaic fragments with no inverse repeats: {}".format((len(noRepeats)/denom)*100))
print("And at least one transposase: {}".format((trans_count(noRepeats)/denom)*100))

import pandas as pd
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Create descriptive figues of incompatability mosaic network")

parser.add_argument("Master", help="Master table file name")
parser.add_argument("Inc_Link", help="Name of links file with incompatability information")
parser.add_argument("Output", help="Desired name of output file")

args = parser.parse_args()

mdf = pd.read_table(args.Master, sep="\t", index_col=0)
incGroup = mdf.groupby(['Best_rep_name']).count()
commonInc = incGroup[incGroup['locus'] >= 20].index
bigIncGroup = mdf[mdf['Best_rep_name'].isin(commonInc)].groupby(['Best_rep_name', 'character']).count()

outframe = bigIncGroup['locus'].unstack().fillna(0.0)
outframe.rename(columns={'static': 'non-mosaic', "character": ""}, inplace=True)
del outframe.columns.name

figs, axes = plt.subplots(nrows=2, ncols=1)

outframe.plot(kind='bar', stacked=True, color=['c', 'm'], ax=axes[0])

columnNames = ['qlocus', 'qstart', 'qend', 'slocus', 'sstart', 'send', 'lineColor']
ldf = pd.read_table(args.Inc_Link, sep="\t", names=columnNames)
interInc = ldf[ldf['lineColor'] == 'color=red']['slocus'].tolist()

mdf['inter_inc'] = mdf['locus'].map(lambda x: (x in interInc))
interIncGroup = mdf[mdf['Best_rep_name'].isin(commonInc)].groupby(['Best_rep_name', 'inter_inc']).count()
outframe2 = interIncGroup['locus'].unstack().fillna(0.0)
outframe2.rename(columns={True: "Exchanged with other incompatibility group", "inter_inc": "",
                          False: "No exchange outside incompatability group"}, inplace=True)
del outframe2.columns.name

outframe2[["Exchanged with other incompatibility group", "No exchange outside incompatability group"]].plot(kind='bar', stacked=True, color=['c', 'm'], ax=axes[1])

#plt.xlabel("Plasmid Inc Type")
axes[0].set(xlabel="Plasmid Inc Type", ylabel="Number of Plasmids")
axes[1].set(xlabel="Plasmid Inc Type", ylabel="Number of Plasmids")
#plt.ylabel("Number of Plasmids")

plt.tight_layout()
plt.savefig(args.Output)

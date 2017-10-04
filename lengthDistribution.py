import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_table("MasterTable.txt", sep="\t", index_col=0)

binBoundaries = range(0, 2200000, 50000)


df.replace(to_replace={'character': 'static'}, value={'character': 'non-mosaic'}, inplace=True)
fig, ax = plt.subplots(nrows=1, ncols=1)
xcut, xbins = pd.cut(df['length'], bins=binBoundaries, retbins=True, right=False, precision=2)
binSeries = pd.Series(xbins[1:]).map(float).round(2)
xcut = pd.cut(df['length'], bins=binBoundaries, right=False, labels=binSeries)
allBins = binSeries.rename('length').to_frame()
x = df.groupby(['character', xcut]).size().unstack('character', fill_value=0)
y = pd.merge(x, allBins, right_on='length', left_index=True, how='right').fillna(0)
del y['length']
y.plot.bar(ax=ax, color=['c', 'm'])
ax.set_xticklabels(binSeries.map(lambda x: int(x/1000)))
ax.set_ylabel("Number of Plasmids")
ax.set_xlabel("Plasmid size (kbp)")
plt.gca().set_yscale("log")
plt.tight_layout()
plt.savefig("LengthDistribution.png")

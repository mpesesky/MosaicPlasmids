import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import numpy as np

par = ap.ArgumentParser(description="Create a histogram from a single column from the input table")

par.add_argument("Table", help="Filename of the input table")
par.add_argument("Column", help="Name of the column to make the histogram out of")
par.add_argument("-b", "--numbins", metavar="X", type=int, default=20, help="Split into X even bins (default)")
par.add_argument("-M", "--binMax", metavar="M", type=float, default=0, help="Explicitly set upper end of bin range")
par.add_argument("-m", "--binMin", metavar="m", type=float, default=0, help="Explicitly set lower end of bin range")
par.add_argument("-s", "--binSize", metavar="S", type=float, default=10, help="Set bin size (must also set max)")
par.add_argument("-g", "--groupby", type=str, default=None, help="Set a groupby column")
par.add_argument("-l", "--label", type=str, default=None, help="x-axis label")
par.add_argument("-o", "--output", type=str, default=None, help="Output filename")
par.add_argument("-L", "--logscale", action='store_true', help="Y-axis in log scale")
par.add_argument("-c", "--numCutoff", metavar="L", type=int, default=0, help="Exclude plasmids < L genes")

args = par.parse_args()

df = pd.read_table(args.Table, sep="\t", index_col=0)
df.replace(to_replace={args.groupby, 'static'}, value={args.groupby, 'non-mosaic'}, inplace=True)

if args.numCutoff > 0:
    df = df[df['CDS'] >= args.numCutoff]

df = df[df['length'] >= 500]

if args.binMax == 0:
    binBoundaries = args.numbins
else:
    binBoundaries = np.arange(args.binMin, args.binMax, args.binSize)

fig, ax = plt.subplots(nrows=1, ncols=1)

if args.groupby is not None:
    xcut, xbins = pd.cut(df[args.Column], bins=binBoundaries, retbins=True, right=False, precision=2)
    binSeries = pd.Series(xbins[1:]).map(float).round(2)
    xcut = pd.cut(df[args.Column], bins=binBoundaries, right=False, labels=binSeries)
    allBins = binSeries.rename(args.Column).to_frame()
    x = df.groupby([args.groupby, xcut]).size().unstack(args.groupby, fill_value=0)
    y = pd.merge(x, allBins, right_on=args.Column, left_index=True, how='right').fillna(0)
    del y[args.Column]
    y.plot.bar(ax=ax, color=['c', 'm'])
    ax.set_xticklabels(binSeries)
else:
    x = df[args.Column]
    x.plot.hist(alpha=0.5, bins=binBoundaries, legend=True, ax=ax, stacked=False)

ax.set_xlabel(args.label)
ax.set_ylabel("Plasmid Count")
plt.tight_layout()

if args.logscale:
    plt.gca().set_yscale("log")

if args.output is not None:
    plt.savefig(args.output)
else:
    plt.savefig("{}_hist.png".format(args.Column))

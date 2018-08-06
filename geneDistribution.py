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
par.add_argument("-C", "--color", type=str, default='cyan', help="Matplotlib color code for plot")
par.add_argument("-p", "--p_val", type=str, default="", help="P-value to print on chart for comparison")

args = par.parse_args()

df = pd.read_table(args.Table, sep="\t", index_col=0)


if args.numCutoff > 0:
    df = df[df['CDS'] >= args.numCutoff]

df = df[df['length'] >= 500]
df.replace({'character':{'static':'non-mosaic'}}, inplace=True)

if args.binMax == 0:
    binBoundaries = args.numbins
else:
    binBoundaries = np.arange(args.binMin, args.binMax, args.binSize)

if args.groupby is not None:
    groupList = list(set(df[args.groupby].tolist()))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 6))
    data = []
    dataMax = df[args.Column].max()
    for x in range(len(groupList)):
        data.append(df[args.Column][df[args.groupby] == groupList[x]])
    plt.boxplot(data, labels=groupList, boxprops=dict(color=args.color),
                flierprops=dict(marker='.', markerfacecolor='white', markeredgecolor='black'),
                whiskerprops=dict(linestyle='-', color=args.color), medianprops=dict(color='black'))
    y, h, col = dataMax + 0.025, 0.025, 'k'
    plt.plot([1, 1, 2, 2], [y, y+h, y+h, y], lw=1.5, c=col)
    pValTxt = "p-value {}".format(args.p_val)
    plt.text(1.5, y+h, pValTxt, ha='center', va='bottom', color=col)
    ax.set_ylabel(args.label)
    plt.ylim(ymax=dataMax+0.1)

else:
    fig, ax = plt.subplots(nrows=1, ncols=1)
    x = df[args.Column]
    x.plot.hist(alpha=0.5, bins=binBoundaries, legend=True, ax=ax, stacked=False)
    ax.set_xlabel(args.label)
    ax.set_ylabel("Plasmid Count")
    if args.logscale:
        plt.gca().set_yscale("log")

plt.tight_layout()



if args.output is not None:
    plt.savefig(args.output)
else:
    plt.savefig("{}_hist.png".format(args.Column))

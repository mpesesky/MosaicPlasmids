import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Get fraction of mosaic plasmids by species")

parser.add_argument("TableFile", help="Data table file")
parser.add_argument("-c", "--cutoff", default=10, type=int, help="Minimum number of plasmids per species in results")
parser.add_argument("-o", "--output", default="MosaicTaxonomy.tsv", help="Name of output tsv file (default = "
                                                                         "MosaicTaxonomy.tsv")
parser.add_argument("-p", "--outplt", default="MosaicTaxonomy.png", help="Name of output png file")
parser.add_argument("-g", "--intergenus", default=None, help="Name of intergenus file for coloring")

args = parser.parse_args()

df = pd.read_csv(args.TableFile, sep="\t")
df = df[df['length'] >= 500]

taxonDF = pd.DataFrame({'count': df.groupby(['Genus', 'character']).size()}).reset_index()

taxonPivot = taxonDF.pivot(index='Genus', columns='character', values='count')

taxonPivot.fillna(0, inplace=True)

taxonBaseline = taxonPivot[(taxonPivot.mosaic + taxonPivot.static) >= args.cutoff]

taxonBaseline['Total Plasmid Sequences'] = taxonBaseline.mosaic + taxonBaseline.static
taxonBaseline['Fraction mosaic'] = taxonBaseline.mosaic.map(float)/taxonBaseline['Total Plasmid Sequences'].map(float)

taxonBaseline.to_csv(args.output, sep="\t")

if args.intergenus is not None:
    taxonBaseline['Genus'] = taxonBaseline.index
    intergenus = pd.read_table(args.intergenus, sep="\t")
    coloreddf = pd.merge(taxonBaseline, intergenus, how='left', left_on='Genus', right_on='organisms').fillna(0)
    outplt = coloreddf.plot.scatter(x='Total Plasmid Sequences', y='Fraction mosaic', marker='o', color=coloreddf['proportion_intergenus_plasmids'], cmap=plt.get_cmap('cool'))
    outplt.set_xlim(xmin=-10)
    outplt.set_ylim(ymin=-0.01)
    fig = plt.gcf()
    cbar = fig.get_axes()[1]
    cbar.set_ylabel("Intergenus fraction of mosaic plasmids")
else:
    outplt = taxonBaseline.plot.scatter(x='Total Plasmid Sequences', y='Fraction mosaic', marker='o', color='c')

plt.tight_layout()
plt.savefig(args.outplt)

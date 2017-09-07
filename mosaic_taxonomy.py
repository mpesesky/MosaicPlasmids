import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Get fraction of mosaic plasmids by species")

parser.add_argument("TableFile", help="Data table file")
parser.add_argument("-c", "--cutoff", default=10, type=int, help="Minimum number of plasmids per species in results")
parser.add_argument("-o", "--output", default="MosaicTaxonomy.tsv", help="Name of output tsv file (default = "
                                                                         "MosaicTaxonomy.tsv")

args = parser.parse_args()

df = pd.read_csv(args.TableFile, sep="\t")

taxonDF = pd.DataFrame({'count': df.groupby(['Genus', 'character']).size()}).reset_index()

taxonPivot = taxonDF.pivot(index='Genus', columns='character', values='count')

taxonPivot.fillna(0, inplace=True)

taxonBaseline = taxonPivot[(taxonPivot.mosaic + taxonPivot.static) >= args.cutoff]

taxonBaseline.to_csv(args.output, sep="\t")

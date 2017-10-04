import pandas as pd
import scipy.stats as spy
import matplotlib.pyplot as plt
import math

def import_table(tableFile):
    df = pd.read_table(tableFile, sep="\t", index_col=0)
    return df

def test_gene_presence(plasmidTable, geneTable):
    gene_presence = plasmidTable.isin(geneTable['plasmid ID'].tolist())
    plasmidTF = gene_presence['locus']
    plasmidTF.name = "gene"
    outDF = pd.concat([plasmidTable, plasmidTF], axis=1)
    return outDF

def test_gene_counts(plasmidTable, geneTable):
    geneCounts = geneTable['plasmid ID'].value_counts().to_dict()
    plasmidTable['gene counts'] = plasmidTable['locus'].map(geneCounts).map(float)
    plasmidTable.fillna(value={'gene counts':0}, inplace=True)
    lengthWeights = plasmidTable['length']/200.0
    plasmidTable['GTHP'] = plasmidTable['gene counts']/lengthWeights
    return plasmidTable

def test_gene_proportions(plasmidTable, geneTable):
    geneCounts = geneTable['plasmid ID'].value_counts().to_dict()
    plasmidTable['gene counts'] = plasmidTable['locus'].map(geneCounts).map(float)
    plasmidTable.fillna(value={'gene counts': 0}, inplace=True)
    plasmidTable["Gene_proportion"] = plasmidTable['gene counts']/plasmidTable['CDS']
    return plasmidTable

def gene_count_by_column(plasmidGeneDF, column):
    outDF = plasmidGeneDF.groupby(column, as_index=False)['GTHP'].mean()
    stdS = plasmidGeneDF.groupby(column)['GTHP'].std()
    stdDict = stdS.rename("STDEV").to_dict()
    outDF['STDEV'] = outDF[column].map(stdDict)
    return outDF


def gene_count_bins(plasmidGeneDF, column, binNum):
    bins = pd.cut(plasmidGeneDF['length'], binNum, labels=False)
    plasmidGeneDF['Bin'] = pd.Series(bins, index=plasmidGeneDF.index)
    meanDF = plasmidGeneDF.groupby([column, 'Bin'], as_index=False)['gene counts'].mean()
    stdS = plasmidGeneDF.groupby([column, 'Bin'], as_index=False)['gene counts'].describe()
    stdDF = stdS['std'].to_frame().reset_index(level=[column, 'Bin'])
    print(stdDF)
    outDF = meanDF.merge(stdDF, on=[column, 'Bin'])
    return outDF


def gene_by_column(plasmidGeneDF, column):
    outDF = pd.DataFrame()
    countSeries = plasmidGeneDF[plasmidGeneDF['gene'] == True][column].value_counts()
    countSeries.name = "Present"
    outDF = outDF.append(countSeries)
    countSeries = plasmidGeneDF[plasmidGeneDF['gene'] == False][column].value_counts()
    countSeries.name = "Absent"
    outDF = outDF.append(countSeries)
    return outDF.fillna(0).transpose()

def gene_proportion(plasmidGeneDF, column):
    outDF = plasmidGeneDF.groupby(column, as_index=False)['Gene_proportion'].mean()
    stdS = plasmidGeneDF.groupby(column)['Gene_proportion'].std()
    stdDict = stdS.rename("STDEV").to_dict()
    outDF['STDEV'] = outDF[column].map(stdDict)
    return outDF

def mann_whitneyu(plasmidGeneDF, testCol):
    mosaic = plasmidGeneDF[plasmidGeneDF['character'] == 'mosaic'][testCol]
    static = plasmidGeneDF[plasmidGeneDF['character'] == 'static'][testCol]

    stat, pval = spy.mannwhitneyu(mosaic, static, use_continuity=False)
    return pval


def calc_expected(plasType, df):
    return (float(df.loc[plasType].sum()) * float(df.loc[[('mosaic', True), ('static', True)]].sum()))/float(df.sum())

def bar_plot(df, xColumn, yColumn, outfile):
    groupCounts = df.groupby(xColumn).agg({yColumn: 'value_counts'})
    groupProportion = groupCounts.groupby(level=0).apply(lambda x: x / float(x.sum()))
    outplt = outfile + ".png"
    propDF = groupProportion.loc[[('mosaic', True), ('static', True)]]
    errDF = groupCounts.loc[[('mosaic', True), ('static', True)]]
    propDF = propDF.reset_index(level=1, drop=True).rename({'static': 'non-mosaic'})
    errDF = errDF.reset_index(level=1, drop=True) .rename(index={'static': 'non-mosaic'}, columns={'gene': 'count'})
    plotDF = propDF.join(errDF)
    plotDF['err'] = ((plotDF['gene']*(1-plotDF['gene']))/plotDF['count']).map(math.sqrt)

    plotDF.plot.bar(y='gene', yerr='err', color='c')
    plt.tight_layout()
    plt.savefig(outplt)

    expected = [calc_expected('mosaic', groupCounts), calc_expected('static', groupCounts)]
    print(expected)
    print(errDF)
    print(groupCounts.sum())
    testStat, pval = spy.stats.chisquare(errDF, f_exp=expected)
    print("Chi-sq results: the test statistic is {} and the p-value is {}\n".format(testStat[0], pval[0]))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Make useful charts")

    parser.add_argument("PlasmidTable", help="Table of Plasmid Information")
    parser.add_argument("GeneTable", help="Table of gene information")
    parser.add_argument("-c", "--character", action='store_true', help="Compare by plasmid character")
    parser.add_argument("-g", "--genus", action='store_true', help="Compare by genus")
    parser.add_argument("-p", "--presence", action='store_true', help="Compare gene presence/absence instead of counts")
    parser.add_argument("-o", "--outfile", type=str, default="Gene_by_plasmid.txt", help="Output file name")
    parser.add_argument("-b", "--bins", type=int, metavar="X", default=0, help="Divide count output into X bins")
    parser.add_argument("-P", "--proportion", action='store_true', help="Compare to total gene count per plasmid")
    parser.add_argument("-m", "--mannTest", action='store_true', help="Perform Mann-Whitney U test (character only)")
    parser.add_argument("-t", "--geneThreshold", type=int, default=3, help="Minimum genes per plasmid")

    args = parser.parse_args()

    plasmidDF = import_table(args.PlasmidTable)
    plasmidDF = plasmidDF[plasmidDF['length'] >= 500]
    plasmidDF = plasmidDF[plasmidDF['CDS'] >= args.geneThreshold]
    geneDF = import_table(args.GeneTable)

    if args.presence:
        plasmidDF = test_gene_presence(plasmidDF, geneDF)
    elif args.proportion:
        plasmidDF = test_gene_proportions(plasmidDF, geneDF)
    else:
        plasmidDF = test_gene_counts(plasmidDF, geneDF)

    plasmidDF.to_csv(args.outfile, sep="\t")

    if args.mannTest == True:
        print(mann_whitneyu(plasmidDF, "gene counts"))

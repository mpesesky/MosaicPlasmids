import pandas as pd

links = pd.read_table("links.txt", sep="\t", header=None, names=['qlocus', 'qstart', 'qend', 'slocus', 'sstart', 'send', 'link'], index_col=False)
organisms = pd.read_table("organisms.txt", sep="\t", header=None, index_col=0).to_dict()[1]
orgLoci = organisms.keys()
for locus in orgLoci:
    genus = organisms[locus].split(" ")[0]
    organisms[locus] = genus
#print("Organism dictionary check: {}".format(organisms['NC_010996']))

loci = pd.DataFrame(pd.concat([links['qlocus'], links['slocus']]).unique(), columns=['locus'])
#print(loci)
loci['organisms'] = loci['locus'].map(organisms)
#print("Loci check: {}".format(loci[loci['locus'] == 'NC_010996']))
output = pd.DataFrame(loci.groupby('organisms').locus.nunique())
output['organisms'] = output.index
#print("Output check 1: {}".format(output[output['organisms'] == 'Rhizobium']))
links['qGenus'] = links['qlocus'].map(organisms)
links['sGenus'] = links['slocus'].map(organisms)

genera = loci['organisms'].unique().tolist()
genusDict = {}
interDict = {}
for org in organisms.keys():
    genusDict[organisms[org]] = 0
    interDict[organisms[org]] = 0

print("Adding links per genus")
for index, row in links.iterrows():
    if (row['qGenus'] == row['sGenus']):
        genusDict[row['qGenus']] += 1
    else:
        genusDict[row['qGenus']] += 1
        genusDict[row['sGenus']] += 1
        interDict[row['sGenus']] += 1
        interDict[row['qGenus']] += 1
print("Final calculations")

interGenusLinks = links[links['qGenus'] != links['sGenus']]
interLoci = pd.DataFrame(pd.concat([interGenusLinks['qlocus'], interGenusLinks['slocus']]).unique(), columns=['locus'])
interLoci['organisms'] = interLoci['locus'].map(organisms)
interPlasmidCounts = pd.DataFrame(interLoci.groupby('organisms').locus.nunique())
interPlasmidCounts['organisms'] = interPlasmidCounts.index
interPlasmidCounts.rename(columns={'locus': 'intergenus_plasmids'}, inplace=True)

output['total_links'] = output['organisms'].map(genusDict)
output['intergenus_links'] = output['organisms'].map(interDict)
output['average_links'] = output['total_links']/output['locus']
output['proportion_intergenus'] = output['intergenus_links']/output['total_links']
output2 = output.merge(interPlasmidCounts, on='organisms', how='outer').fillna(0)
output2['proportion_intergenus_plasmids'] = output2['intergenus_plasmids'].map(float)/output2['locus'].map(float)

output2.to_csv("Link_stats.tsv", sep='\t')

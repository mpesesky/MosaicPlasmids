import argparse

def get_bioSamples(gbFile):
    from Bio import SeqIO
    sampleDict = {}
    for record in SeqIO.parse(gbFile, "genbank"):
        locus = record.id
        dblinks = record.dbxrefs
        if len(dblinks) > 1:
            for db in dblinks:
                fields = db.split(":")
                if fields[0] == "BioSample":
                    sampleDict[locus] = fields[1]
    return sampleDict

def get_environments(sampleDict):
    from Bio import Entrez
    import xml.etree.ElementTree as ET
    Entrez.email = "mwpesesky@gmail.com"
    envDict = {}
    for locus in sampleDict.keys():
        biosampleCall = Entrez.efetch(db='biosample', id=sampleDict[locus], retmode='xml')
        root = ET.fromstring(biosampleCall.read())
        biosampleCall.close()
        found = False
        try:
            attribs = root[0].findall('Attributes')
        except IndexError:
            print("XML Error: {}".format(locus))
            continue
        for attrib in attribs[0]:
            if (attrib.get('attribute_name') == 'Isolation Site') or (attrib.get('attribute_name' == 'isolation site')):
                envDict[locus] = attrib.text
                found = True
                continue
        for attrib in attribs[0]:
            if (attrib.get('attribute_name') == 'Isolation Source') or (attrib.get('attribute_name') == 'isolation source'):
                envDict[locus] = attrib.text
                found = True
                continue
        for attrib in attribs[0]:
            if ('environment' in attrib.get('attribute_name')) or ('Environment' in attrib.get('attribute_name')):
                envDict[locus] = attrib.text
                found = True
                continue

        if not found:
            print("No env: {}".format(locus))
    return envDict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the Environmental sources from a GBfile")

    parser.add_argument("GenBankFile", help="File to get samples from")
    parser.add_argument("-o", "--output", type=str, default="Environments.txt", help="Output file name ("
                                                                                     "default=Environments.txt")

    args = parser.parse_args()

    bioSamples = get_bioSamples(args.GenBankFile)

    environments = get_environments(bioSamples)

    outfile = open(args.output, 'w')

    for locus in environments.keys():
        outfile.write("{}\t{}\t{}\n".format(locus, environments[locus], bioSamples[locus]))
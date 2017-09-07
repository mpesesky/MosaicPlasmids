#! /usr/bin/env python

#Adapted from Perl by Mitch Pesesky (original by Alejandro Reyes)

# Receives a Fasta File and print all the sequences in only one line each
#
# Usage: Fasta_one_line.py "Fasta File" > "output"
#
# Input: Fasta File
# Output: Modified Fasta File
# Created: Aug 11 / 08.
# Last-updated: Feb 01 / 16.

def one_line_o(infile):

    f = open(infile)

    seq = ""

    for line in f:
        if line[0] == '>':
            print(seq.rstrip("\n"))
            seq = ""
            print(line.rstrip("\n"))
        else:
            seq = seq.rstrip("\n") + line

    print(seq.rstrip("\n"))
    f.close()
    return True

def one_line_d(infile):
    
    f = open(infile)

    seq = ""
    output = {}
    header = "X"

    for line in f:
        if line[0] == '>' and header == "X":
            header = line.rstrip()
        elif line[0] == '>':
            output[header] = seq
            header = line.rstrip()
            seq = ""
        else:
            seq = seq + line.rstrip()
    if seq != "":
        output[header] = seq
    f.close()
    return output

def blast_dict(infile):
    f = open(infile)

    seq = ""
    output = {}
    header = "X"

    for line in f:
        if line[0] == '>' and header == "X":
            header = line.lstrip(">").split(" ")[0]
        elif line[0] == '>':
            output[header] = seq
            header = line.lstrip(">").split(" ")[0]
            seq = ""
        else:
            seq = seq + line.rstrip()
    if seq != "":
        output[header] = seq
    f.close()
    return output

def one_line_l(infile):
    
    f = open(infile)

    seq = ""
    output = []
    first = True

    for line in f:
        if line[0] == '>' and first:
            output.append(line.rstrip("\n"))
            first = False
        elif line[0] == '>':
            output.append(seq)
            output.append(line.rstrip("\n"))
            seq = ""
        else:
            seq = seq + line.rstrip("\n")

    output.append(seq)
    f.close()
    return output


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        sys.exit("\nUsage: Fasta_one_line.py <input.fasta> > output.fasta\n\n")

    if one_line_o(sys.argv[1]):
        print("Complete")
    else:
        sys.exit("Error, could not complete request. Check script and input file")

            





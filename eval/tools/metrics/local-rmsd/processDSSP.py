#!/usr/bin/python

"""
DESCR: Returns the strings of the sequences of amino acids and secondary structures of a PDB file.
       It processes the output from the 'dsspcmbi' program.
AUTHOR: Luis Garreta
DATE:	May, 13, 2010

NOTES:
	It ASSUMES that always both amino and the structue will be in the positions 13, 16 respectivelly.
	This assumption because it was impossible to split the line by a pattern as there are many empty values
"""

import sys

USAGE  = "Returns the strings of the sequences of amino acids and secondary structures of a PDB file.\n"
USAGE += "It processes the output from the 'dsspcmbi' program.\n"
USAGE += "USAGE: python processDSSP.py <pdb filename> [optional output filename]\n"

###############################################################################
if len (sys.argv) == 2:
	None
elif len (sys.argv) == 3:
	outputFilename = sys.argv [2]
	sys.stdout = open (outputFilename, "w")
else:
	print USAGE
	sys.exit (0)

###############################################################################

dsspFilename = sys.argv [1]
dsspFile = open (dsspFilename)

# Skips header lines until reach the values of aminos and structures and others
while "#" not in dsspFile.next():
	None

# Starts to copy the amino and its corresponding structure. It ASSUMES that "always"
# both amino and the structue will be in the positions 13, 16 respectivelly.
sequenceOfAminoNumbers = ""
sequenceOfAminos = ""
sequenceOfStructures = ""

for line in dsspFile:
	aminoNumber = line [9]
	sequenceOfAminoNumbers += aminoNumber

	amino = line [13]
	sequenceOfAminos += amino

	structure = line [16]
	sequenceOfStructures += structure

print sequenceOfAminoNumbers
print sequenceOfAminos
print sequenceOfStructures
	

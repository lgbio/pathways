#!/usr/bin/python

import os, sys
###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
	baseName = os.path.basename (namefile)
	newNamefile = baseName.split(".") [0]

	return newNamefile

###############################################################################
def calculateDssp (inputPDB, outputDssp):
	
	command = "dsspcmbi " + inputPdb + " " + outputDssp
	os.system (command)

###############################################################################
def formatDsspToAminoStructure (inputDssp, outputAminoStructure):

	dsspFile = open (inputDssp)
	# Skips header lines until reach the values of aminos and structures and others
	while "#" not in dsspFile.next():
		None

	# Starts to copy the amino and its corresponding structure. It ASSUMES that "always"
	# both amino and the structue will be in the positions 13, 16 respectivelly.
	sequenceOfAminos = ""
	sequenceOfStructures = ""

	for line in dsspFile:
		amino = line [13]
		sequenceOfAminos += amino

		structure = line [16]
		sequenceOfStructures += structure

	
	outputFastaFile = open (outputAminoStructure, "w")
	outputFastaFile.write (sequenceOfAminos + "\n")
	outputFastaFile.write (sequenceOfStructures)
	outputFastaFile.close()

###############################################################################
def calculateRateOfSecondaryStructures (inputAminoStructure):
	aminoStructureFile = open (inputAminoStructure)
	aminos = aminoStructureFile.next()
	structures = aminoStructureFile.next()

	sizeStructures = len (structures)

	countOfStructures = 0
	countOfCoils = 0
	for str in structures:
		if str == " ":
			countOfCoils += 1
		else: #if s == "H" or s == "B" or s == "E" or s == "G" or str == "I":
			countOfStructures += 1

	return float (countOfStructures) / sizeStructures

###############################################################################
def run (inputPdb):
	dssp = name (inputPdb) + ".dssp"
	print ">>>", dssp
	calculateDssp (inputPdb, dssp)

	aminoStructure = name (dssp) + ".aass"
	formatDsspToAminoStructure (dssp, aminoStructure)

	rateSecondaryStructures = calculateRateOfSecondaryStructures (aminoStructure)

	print rateSecondaryStructures

###############################################################################

if len (sys.argv) == 2:
	inputPdb = sys.argv [1]
	run (inputPdb)
else:
	print USAGE
	sys.exit (0)


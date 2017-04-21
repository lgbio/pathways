#!/usr/bin/python

USAGE = "\
Returns the radius of gyration for a protein using the gromacs g_gyrate (opt. 1)\n\
USAGE: radiusg-lg.py <PDB structure>\n"

import os, sys

stderr = os.getenv ("EVAL_STDERR")
sys.stderr = open (stderr, "a")

###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> RdGy: "
	if (type (message) != list): 
		string +=  str (message)
	else: 
		for i in message: string += str (i) + " "
	
	sys.stderr.write (string+"\n")

###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
	baseName = os.path.basename (namefile)
	newNamefile = baseName.split(".") [0]

	return newNamefile

###############################################################################
# Remove temporal filenames used for calculations
###############################################################################
def clean (listOfTmpFiles):
	for tmpFile in listOfTmpFiles:
		os.system ("rm " + tmpFile)

##################################################################
# RUN
##################################################################
def run (pdbFilename):
	PDBStructure = pdbFilename
	logFilename = name (pdbFilename) + "-rg-tmp.log"
	xvgFilename = name (pdbFilename) + "-rg-tmp.xvg"

	command = "echo 1 | g_gyrate -f %s -s %s -o %s >& %s " % (PDBStructure, PDBStructure, xvgFilename, logFilename)
	log (command)
	os.system (command)

	resultsOfRmsd = open (xvgFilename).readlines()
	radiusOfGyration = -1

	for line in resultsOfRmsd:
		if "#" in line[0] or "@" in line[0]:
			continue
		else:
			radiusOfGyration = float (line.split()[1].strip()) * 10
			break

	clean ([xvgFilename, logFilename])
	return radiusOfGyration

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":
	if len (sys.argv) < 2:
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv [1]

	print run (pdbFilename)


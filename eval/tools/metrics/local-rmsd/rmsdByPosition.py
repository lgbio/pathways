#!/usr/bin/python

"""
GOAL: Calculate the RMSD by position (start, end) of two PDBs
INPUT: Reference PDB, Query PDB, start position, end postion
OUPUT: RMSD Value
"""
USAGE  = "Calculate the RMSD by position (start, end) of two PDBs\n"
USAGE += "USAGE: rmsdByPosition.py <refPdb> <qryPdb> <start> <end>"

import sys, os, subprocess

###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
	return value

###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
	baseName = os.path.basename (namefile)
	newNamefile = baseName.split(".") [0]

	return newNamefile

###############################################################################
def renumber (refPdb, renumberPdb):
	cmm = "convpdb.pl -renumber %s %s > %s" % ("1", refPdb, renumberPdb)
	print cmm
	os.system (cmm)
	
###############################################################################
def cut (inPdb, start, end, outPdb):
	cmm = "convpdb.pl -sel %s:%s %s > %s" % (start, end, inPdb, outPdb)
	print cmm
	os.system (cmm)
	
###############################################################################
def rmsd (refPdb, qryPdb):
	print "rmsd " + refPdb + " " + qryPdb
	params = ["rmsd", refPdb, qryPdb]
	value = runProgram (params, ".")

	return value

###############################################################################
def clean (listOfTmpFiles):
	for i in listOfTmpFiles:
		os.system ("rm " + i)
###############################################################################
def rmsdByPosition (refPdb, qryPdb, start, end):
	renumberPdb = name (refPdb) + "-renumber.pdb"
	renumber (refPdb, renumberPdb)

	cutRef = name (renumberPdb) + "-cut.pdb"
	cut (renumberPdb, start, end, cutRef)

	cutQuery = name (qryPdb) + "-cut.pdb"
	cut (qryPdb, start, end, cutQuery)

	value = rmsd (cutRef, cutQuery)

	#clean ([renumberPdb, cutRef, cutrenumber, "*.rmsd"])

	return value
###############################################################################
# Main
###############################################################################
if __name__ == "__main__":

	if len (sys.argv) != 5:
		print USAGE
		sys.exit (0)
	else:
		refPdb = sys.argv [1]
		qryPdb = sys.argv [2]
		startPos = sys.argv [3]
		endPos = sys.argv [4]

		value = rmsdByPosition (refPdb, qryPdb, startPos, endPos)

		print ">>>", value

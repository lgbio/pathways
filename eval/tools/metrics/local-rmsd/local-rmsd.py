#! /usr/bin/python 

import os, sys, subprocess
import uuid

"""
LOG:
	2016/12/16:  Added handle of stderr messages
"""

USAGE = "\
GOAL:   Calculates the local RMSD of a reference PDB vs. query PDB (Average of SS).\n\
INPUT:  Reference PDB and query PDB.\n\
OUTPUT: Average RMSD of secondary structures. \n\
USAGE:  local_rmsd.py <reference PDB> <query PDB> \n"

TMPLABEL = "tmp_LR_%s_" % str (uuid.uuid4())

#-------------------------------------------------------------
# Main 
#-------------------------------------------------------------
def main (args):
	global TMPLABEL
	if len (args) != 3:
		print USAGE
		sys.exit (0)
	else:
		targetPDB = args [1]
		referencePDB = args [2]
		TMPLABEL += extractName (referencePDB) + "_"
		log ("%s, %s" % (targetPDB, referencePDB))
		localRmsdValue = localRmsd (targetPDB, referencePDB)
		print localRmsdValue

#-------------------------------------------------------------
# Calculate the local RMSD
#-------------------------------------------------------------
def localRmsd (targetPDB, refPDB):
	matchPdb = matchResidueNumbering (targetPDB, refPDB)
	refInfo = getAminoStructureInfo (refPDB)
	refStructures = getLocalStructures (refInfo["aminoStructure"], refInfo["startResidueNumber"])
	localRmsdValue = calculateAverageRmsd (targetPDB, matchPdb, refStructures)

	clean ([matchPdb])
	return localRmsdValue


###############################################################################
# Call external rmsd 
###############################################################################
def rmsd (targetPDB, refPDB, start, end):
	params = ["rmsd", targetPDB, refPDB, str(start), str (end)]
	value = runProgram (params, ".")
	return float (value)

###############################################################################
# Call external dssp 
###############################################################################
def calculateDssp (PDB):
	dsspName = TMPLABEL + extractName (PDB) + ".dssp"
	cmm = "dssp %s %s" % (PDB, dsspName)
	log (cmm)
	params = ["dssp", PDB, dsspName]
	value = runProgram (params, ".")

	return dsspName
###############################################################################
# Calculate the DSSP of "PDB" structure and get from it information:
# sequence A+B, start residue number, end residue number and size of the sequence
###############################################################################
def getAminoStructureInfo (PDB):
	dsspName = calculateDssp (PDB)
	dsspFile = open (dsspName)
	# Skips header lines until reach the values of aminos and structures and others
	while "#" not in dsspFile.next():
		None

	# Starts to copy the amino and its corresponding structure. It ASSUMES that "always"
	# both amino and the structue will be in the positions 13, 16 respectivelly.
	sequenceOfAlphaBeta = ""
	start = ""
	for line in dsspFile:
		if start == "":
			start = line.split ()[1]

		structure = line [16]

		if structure in ['H','G','I']:
			sequenceOfAlphaBeta += 'H'
		elif structure in ['B','E']:
			sequenceOfAlphaBeta += 'B'
		else:
			sequenceOfAlphaBeta += 'C'
	
	end = line.split ()[1]

	info = {}
	info ["aminoStructure"] = sequenceOfAlphaBeta
	info ["startResidueNumber"] = int (start)
	info ["endResidueNumber"] = int (end)
	info ["ziseOfSequence"] = len (sequenceOfAlphaBeta)

	clean ([dsspName])
	return info

###############################################################################
# Get a substructure from a sequence "seq" that starts in the position "pos"
###############################################################################
def getSubStructure (struct, pos, seq, lstStructs, startResidueNumber):
	start = pos 
	pos += 1
	nextStruct = seq [pos]
	size = len (seq)
	while pos < size  and nextStruct == struct: 
		pos += 1
		nextStruct = seq [pos]

	end = pos + 1
	pos -= 1
	lstStructs.append ((start+startResidueNumber-1, end+startResidueNumber-1))

	return pos

###############################################################################
# Get local structures from "aminoStructure" sequence starting in a "startResidueNumber"
# When it detects a starting SS, it calls the getSubStructure 
###############################################################################
def getLocalStructures (aminoStructure, startResidueNumber):
	position=0
	structure=""
	listOfStructures = []
	sizeStructure = len (aminoStructure)
	while position < sizeStructure:
		structure = aminoStructure [position]

		if structure in ["H", "B"]:
			position = getSubStructure (structure, position, aminoStructure, listOfStructures, startResidueNumber)

		position += 1

	return listOfStructures

###############################################################################
# Match the residue numbering of a reference PDB to a query PDB
###############################################################################
def matchResidueNumbering (targetPDB, refPDB):
	matchPdb = "%s_%s_%s_%s" % (TMPLABEL,extractName(targetPDB),extractName(refPDB),"_match.pdb")
	cmm = "convpdb.pl -match %s %s > %s" % (targetPDB, refPDB, matchPdb)
	os.system (cmm)

	return matchPdb
	
###############################################################################
# 
###############################################################################
def calculateAverageRmsd (matchPdb, refPdb, listOfStructures):
	listOfRmsds =[]

	sumRmsd = 0
	sumSizes = 0
	for s in listOfStructures:
		start = int (s[0])
		end = int (s[1])
		sizeStructure = end - start + 1

		rmsdValue = rmsd (matchPdb, refPdb, start, end) * sizeStructure

		sumRmsd += rmsdValue
		sumSizes += sizeStructure

	if sumSizes == 0:
		totalRmsd = 99
	else:
		totalRmsd = sumRmsd / float (sumSizes)

	return totalRmsd

###############################################################################
# Remove temporal filenames used for calculations
###############################################################################
def clean (listOfTmpFiles):
	for tmpFile in listOfTmpFiles:
		os.system ("rm " + tmpFile)
	
###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> LRMSD: "
	if type (message) != list:
		string +=  str (message)
	else:
		for i in message: string += str (i) + " "
	
	sys.stderr.write (string+"\n")

#------------------------------------------------------------
# Extract only the name without path and extension
#------------------------------------------------------------
def extractName (fullName):
	return os.path.splitext (os.path.basename (fullName))[0]

#-------------------------------------------------------------
# Define the  output for errors and log messages
#-------------------------------------------------------------
def defineMessagesOutput ():
	stderr = os.getenv ("EVAL_STDERR")
	if stderr == None:
		sys.stderr = sys.stdout
	else:
		sys.stderr = open (stderr, "a")

###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE, stderr=sys.stderr).communicate()[0]
	return value

###############################################################################
# Main
###############################################################################
if __name__ == "__main__":
	defineMessagesOutput ()
	main (sys.argv)



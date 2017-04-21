#!/usr/bin/python

import os, sys
import tarfile       # For extractTrajectory function
import subprocess    # For runProgram function
from multiprocessing import Pool # For MultiProcessing Class For Multiprocessing

"""
Evaluate diferent metrics (as scripts) to a single pathway or
a directory of pathways
"""
USAGE = "eval01_main.py <-f|-d> <single pathway|dir of pathways> [outut dir] [nCpus]"
#-----------------------------------------------------------
#-----------------------------------------------------------
def main (args):
	if len (args) < 3:
		print ">", USAGE
		sys.exit (0)
		
	if args [1] == "-f":
		fullPathwayName = args [2]
		print "Single path", fullPathwayName
		outputDir = os.getcwd ()
		evalSinglePathway (fullPathwayName, outputDir)
	elif args [1] == "-d":
		inputDirPathways = args [2]
		outputDir   = args [3]
		nCpus       = int (args [4])
		evalDirPathways (inputDirPathways, outputDir, nCpus)

#-----------------------------------------------------------
#-----------------------------------------------------------
def evalDirPathways (inputDirPathways, outputDir, nCpus=1):
	checkExistingDir (outputDir)
	os.mkdir (outputDir)

	inputFiles = sorted ([x for x in os.listdir (inputDirPathways) if ".tgz" in x])
	inputFilesFullPath = ["%s/%s" % (inputDirPathways, x) for x in inputFiles]

	for fullPath in inputFilesFullPath:
		evalSinglePathway (fullPath, outputDir, nCpus)

	"""
	params = zip (inputFilesFullPath, [outputDir]*len (inputFiles))
	pool = Pool (processes=nCpus)
	pool.map (evalSinglePathway_Start, params)
	pool.close ()
	pool.join ()
	"""

#-----------------------------------------------------------
# Evaluates the 11 properties to a single pathway (.tgz)
#-----------------------------------------------------------
def evalSinglePathway_Start (args):
	evalSinglePathway (*args)

def evalSinglePathway (fullPathwayName, outputDir, nCpus=1):
	print ">>> Evaluating properties on %s ..." % fullPathwayName
	headersOutputFile = "PDB NC CO RG HB AS RM LR RC RA DF SS" 

	tmpDir, refPdbFullpath = extractTrajectoryTgz (fullPathwayName, outputDir)

	inputFiles = sorted (os.listdir (tmpDir))
	inputFilesFullPath = ["%s/%s" % (tmpDir, x) for x in inputFiles]
	params = zip (inputFilesFullPath, [refPdbFullpath]*len (inputFiles))

	pool = Pool (processes=nCpus)
	results = pool.map (evalSinglePdbFile_Star, params)

	# Write to file
	resultsFilename = "%s/%s.values" % (outputDir, extractName (fullPathwayName))
	resultsFile = open (resultsFilename, "w")
	resultsFile.write ("\t".join (headersOutputFile.split()) + "\n")

	for row in results:
		resultsFile.write (row + "\n")
	resultsFile.close ()

	#os.system ("rm -rf %s &" % tmpDir)

#-----------------------------------------------------------
#-----------------------------------------------------------
def evalSinglePdbFile_Star (args):
	return evalSinglePdbFile (*args)

def evalSinglePdbFile (pdbFileFullpath, refPdbFullpath):
	valuesList = [extractName (pdbFileFullpath)]

	cmm = 'native-contacts %s %s' % (pdbFileFullpath, refPdbFullpath)
	valuesList.append (executeCommand (cmm))
	cmm = "contact-order %s" % pdbFileFullpath
	valuesList.append (executeCommand (cmm))
	cmm = "radius-gyration %s" % pdbFileFullpath
	valuesList.append (executeCommand (cmm))
	cmm = "hydrogen-bonds %s" % pdbFileFullpath
	valuesList.append (executeCommand (cmm))
	cmm = "accessible-surface %s" % pdbFileFullpath
	valuesList.append (executeCommand (cmm))
	cmm = 'rmsd %s %s' % (pdbFileFullpath, refPdbFullpath)
	valuesList.append (executeCommand (cmm))
	cmm = 'local-rmsd %s %s' % (pdbFileFullpath, refPdbFullpath)
	valuesList.append (executeCommand (cmm))
	cmm = 'secondary-structures-correct %s %s' % (pdbFileFullpath, refPdbFullpath)
	valuesList.append (executeCommand (cmm))
	cmm = 'secondary-structures-any %s' % (pdbFileFullpath)
	valuesList.append (executeCommand (cmm))
	cmm = 'degrees_freedom %s' % (pdbFileFullpath)
	valuesList.append (executeCommand (cmm))
	cmm = 'structural-similarity %s %s' % (pdbFileFullpath, refPdbFullpath)
	valuesList.append (executeCommand (cmm))

	return "\t".join (valuesList)

#-----------------------------------------------------------
# Extract the conformations (Pdbs) of the trajectory and
# get the starting and reference (ending) pdbs.
# NOTE: Assumes that pathway is preprocessed (parasol or antom)
#-----------------------------------------------------------
def extractTrajectoryTgz (fullPathwayName, outputDir):
	pathwayName = extractName (fullPathwayName)
	tmpDir = "%s/tmp_pdbs_%s" % (outputDir, pathwayName)

	tar = tarfile.open (fullPathwayName, "r|gz")
	tar.extractall (path=tmpDir)

	refPdbFullpath = "%s/%s" % (tmpDir, tar.getnames ()[-1])

	return tmpDir, refPdbFullpath

#------------------------------------------------------------
# call a external programm that returns a value running on "workingDir"
#------------------------------------------------------------
def executeCommand (cmm, workingDir="./"):
	print ">>>", cmm,
	try:
		value = subprocess.Popen (cmm.split (" "), cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
		print value,
		return value.strip()
	except:
		log(["runProgram Error:"] + cmm.split())

#-------------------------------------------------------------
# Print a message to the error output stream
#-------------------------------------------------------------
def log (message):
	string=">>> PropEval: "
	if (type (message) != list): 
		string +=  str (message)
	else: 
		for i in message: string += str (i) + " "
	
	sys.stderr.write (string+"\n")

#----------------------------------------------------------------------
# Move dir to old-dir. Used when a directory exists.
#----------------------------------------------------------------------
def checkExistingDir (dir):
	if os.path.lexists (dir):
		headDir, tailDir = os.path.split (dir)
		oldDir = os.path.join (headDir, "old-" + tailDir)
		if os.path.lexists (oldDir):
			checkExistingDir (oldDir)

		os.rename (dir, oldDir)
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


#------------------------------------------------------------
# MAIN
#------------------------------------------------------------

if __name__ == "__main__":
	os.environ ["EVAL_STDERR"] = os.getcwd()+"/errors.log"
	stderr = os.getenv ("EVAL_STDERR")
	#sys.stderr = open (stderr, "w")

	main (sys.argv)


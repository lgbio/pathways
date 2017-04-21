#!/usr/bin/python

"""
 For parasol trajectories:
 	Extract the conformations (Pdbs) of the trajectory and
 	unify their names (e.g. conformation_10 -> conformation_010)
"""

import os, sys
import tarfile       # For extractTrajectory function
from multiprocessing import Pool # For MultiProcessing Class For Multiprocessing

USAGE = sys.argv [0] + " -f|-d <pathway name|dir name>"
#-----------------------------------------------------------
# Main
#-----------------------------------------------------------
def main (args):
	if len (args) < 3:
		print USAGE
		sys.exit (0)

	if args [1] == "-f":
		pathwayName = args [2]
		preprocessSinglePathway (pathwayName)
	elif args [1] == "-d":
		print ">>>", "Processing multiple pathways..."
		pathwaysDir    = args [2]
		outputDir      = args [3]
		nCpus          = int (args [4])
		preprocessMultiplePathways (pathwaysDir, outputDir, nCpus)

#-----------------------------------------------------------
#-----------------------------------------------------------
def preprocessMultiplePathways (pathwaysDir, outputDir, nCpus=1):
	print ">>> pd", pathwaysDir
	inputFiles = [x for x in os.listdir (pathwaysDir) if ".tgz" in x]
	fullInputFiles = ["%s/%s"%(pathwaysDir,x) for x in sorted (inputFiles)]

	pool = Pool (processes=nCpus)
	params = zip (fullInputFiles, [outputDir] * len (fullInputFiles))
	pool.map (preprocessSinglePathway_STAR, params)
	pool.close ()
	pool.join ()

#-----------------------------------------------------------
#-----------------------------------------------------------
def preprocessSinglePathway_STAR (args):
	preprocessSinglePathway (*args)

def preprocessSinglePathway (pathwayName, outputDir="."):
	newPathwayName = "%s/%s_prep.tgz" % (outputDir, extractName (pathwayName))
	tmpDir = "tmp_%s" %  pathwayName.split(".")[0]
	print ">>> pn", newPathwayName, tmpDir

	tar = tarfile.open (pathwayName, "r|gz")
	tar.extractall (path=tmpDir)
	listOfPdbs = filter (lambda x: ".pdb" in x, os.listdir (tmpDir))

	for i, pdb in enumerate (listOfPdbs):
		name = extractName (pathwayName).split("-pdbs")[0]
		number = pdb.split (".")[0].split ("-")[1]
		newName = tmpDir + "/" + name + "-" + number.zfill (3) + ".pdb"

		os.rename ("%s/%s" % (tmpDir, pdb), newName)
		listOfPdbs [i] = newName

	with tarfile.open (newPathwayName, "w|gz") as tar:
		for path in sorted (listOfPdbs):
			print ">>>", path
			tar.add (path,arcname=os.path.basename (path))

	os.system ("rm -rf %s" % tmpDir)

#------------------------------------------------------------
# Extract only the name without path and extension
#------------------------------------------------------------
def extractName (fullName):
	return os.path.splitext (os.path.basename (fullName))[0]

#------------------------------------------------------------
# MAIN
#------------------------------------------------------------
if __name__ == "__main__":
	main (sys.argv)



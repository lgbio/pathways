#!/usr/bin/python

"""
Return the total void volume in a protein. It uses the 'avp' program.
It uses a default probe of 0.5 Angstrom

AUTHOR: Luis Garreta
DATE: Dic 06/2010

"""
USAGE  = "Return different values of rigidity analysis. It uses the 'FIRST' program\n"
USAGE += "USAGE	:  rigidity-analysis.py <Type of analysis flag>  <pdb filename>\n"
USAGE += "INPUT	:  Type of analysis: --RigidClusters|--StressedRegions|--DegreesOfFreedom\n"

import sys, os, uuid

TMPLABEL="tmp_RIGIDITY_%s_" % str (uuid.uuid4())

###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> RIGIDITY: "
	for i in message:
		string += str (i) + " "
	
	sys.stderr.write (string+"\n")

###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
	baseName = os.path.basename (namefile)
	newNamefile = baseName.split(".") [0]

	return newNamefile

##################################################################
# Delete temporal files created during the process
##################################################################
def clean (listOfFiles):
	for file in listOfFiles:
		os.system ("rm " + file)

##################################################################
# Process the output file and get the values for three properties
##################################################################
def getRigidityValues (outFilename):
	with open (outFilename) as outFile:
		#rigidCluster, stressedRegions, degreesOfFreedom = 0, 0, 0
		for line in outFile:
			line = line.lstrip()
			if "Rigid clusters larger" in line:
				rigidCluster = line.split (" ")[0]
			elif "Stressed Regions" in line:
				stressedRegions = line.split (" ")[0]
			elif "degrees of freedom" in line:
				degreesOfFreedom = line.split (" ")[0]
				break

	return rigidCluster, stressedRegions, degreesOfFreedom

###############################################################################
# 
###############################################################################
def run (pdbFilename, rigidityValueType):
	lib = os.getenv ("EVAL_FIRST")
	if lib == None:
		log (["Error loading FIRST LIB"])
		sys.exit (0)

	params = {"LIB":lib,"IN": pdbFilename, "OUT": "/dev/null"}
	command = "FIRST -non -E -1.0 -L %(LIB)s %(IN)s > %(OUT)s" % params
	log ([command])
	os.system (command)

	outFilename = pdbFilename.split (".")[0] + "_results.txt"
	rigidClusters, stressedRegions, degreesOfFreedom = getRigidityValues (outFilename)

	#clean ([outFilename])

	if rigidityValueType == "--RigidClusters":
		return rigidClusters
	elif rigidityValueType == "--StressedRegions":
		return stressedRegions
	elif rigidityValueType == "--DegreesOfFreedom":
		return degreesOfFreedom
##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) != 3: 
		print USAGE
		sys.exit (0)

	rigidityValueType = sys.argv [1]
	pdbFilename = sys.argv[2]

	value = run (pdbFilename, rigidityValueType)
	print value 


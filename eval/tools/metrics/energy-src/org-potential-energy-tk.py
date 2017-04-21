#!/usr/bin/python

"""
AUTHOR: Luis Garreta
DATE: May 07/2010

NOTES:
	It has problems when the protein has two chains. The "pdbxyz" programs requires an
	additional parameter and it will be waiting to get this one.

"""
USAGE  = "Return the 'Total Potential Energy' for a protein using tinker 'analyze' program\n"
USAGE += "USAGE:  potential-energy-tk.py <pdb filename>\n"
USAGE += "INPUT:  the pdb filename\n"
USAGE += "OUTPUT: the value of the potential energy"

import sys, os

# Load enviromental variables (tinker, force field)
TINKER_PATH     = os.getenv ("TINKER_PATH")
PHD_FORCE_FIELD = os.getenv ("PHD_FORCE_FIELD")

if PHD_FORCE_FIELD == "AMBER96":
	FORCE_FIELD = TINKER_PATH + "/params/amber02.prm"
elif PHD_FORCE_FIELD == "OPLSAA":
	FORCE_FIELD = TINKER_PATH + "/params/oplsaa.prm"
else: 
	print "ERROR in potential energy: Invalid Force Field"
	sys.exit()

###############################################################################
# Return the stem name of a relative or full name with some extension
###############################################################################
def name (namefile):
	reverseNamefile = namefile [::-1]
	pos = reverseNamefile.index (".")
	newNamefile = reverseNamefile [pos+1:][::-1]

	return newNamefile

##################################################################
# Delete temporal files created during the process
##################################################################
def clean (listOfFiles):
	for file in listOfFiles:
		os.system ("rm " + file)

###############################################################################
#
###############################################################################
def createLogFilename (pdbFilename):
	logFilename = name (pdbFilename) + "-tmp-energy.log"
	return logFilename

###############################################################################
#
###############################################################################
def createXyzFromPdb (pdbFilename, FORCE_FIELD, logFilename):
	os.system ("pdbxyz " + pdbFilename + " " + FORCE_FIELD + " > " + logFilename)
	xyzFilename = name (pdbFilename) + ".xyz"
	return xyzFilename

###############################################################################
# Create the temporal file including the parameters for 'analyze'program
###############################################################################
def createTinkerParametersFile (xyzFilename, FORCE_FIELD):
	proteinTinkerFilename = name (xyzFilename) + "-tmp-tinker.params"
	proteinFormatFile = open (proteinTinkerFilename, "w")

	proteinFormatFile.write (xyzFilename + "\n")
	proteinFormatFile.write (FORCE_FIELD + "\n")
	proteinFormatFile.write ("E")
	proteinFormatFile.close ()

	return proteinTinkerFilename

##################################################################
# Get the potential energy by preprocessing the output energy file
##################################################################
def getPotentialEnergy (parametersFilename, logFilename):
	os.system ("analyze < " + parametersFilename  + " >> " + logFilename )
	logFile = open (logFilename)
	potentialEnergy = ""
	for line in logFile:
		if "Total Potential Energy :" in line:
			value = line.split ()[4].strip() 
			value = value.replace ("D", "E")
			return value

	return "9999999999999999.999999999999"

##################################################################
# RUN: Create xyz file, parameters file, and preprocess the energy output file
##################################################################
def run (pdbFilename):
	logFilename = createLogFilename (pdbFilename)
	xyzFilename = createXyzFromPdb (pdbFilename, FORCE_FIELD, logFilename )
	parametersFilename = createTinkerParametersFile (xyzFilename, FORCE_FIELD)
	potentialEnergy = getPotentialEnergy (parametersFilename, logFilename)

	#clean ([logFilename, xyzFilename, parametersFilename, name(pdbFilename)+".seq"])
	return potentialEnergy

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv[1]

	print run (pdbFilename)


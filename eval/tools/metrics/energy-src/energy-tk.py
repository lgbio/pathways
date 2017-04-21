#!/usr/bin/python

"""
AUTHOR: Luis Garreta
DATE: May 07/2010

NOTES:
	It has problems when the protein has two chains. The "pdbxyz" programs requires an
	additional parameter and it will be waiting to get this one.

"""
USAGE  = "Return the 'Total Potential Energy' or 'Dipole Moment' for a protein using tinker 'analyze' program\n"
USAGE += "USAGE:  energy <type of analysis> <pdb filename>\n"
USAGE += "INPUT:  type of analysis: --potential|--dipole; and the pdb filename\n"
USAGE += "OUTPUT: the value of the energy analysis"

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
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> ENR: "
	for i in message:
		string += str (i)
	
	sys.stderr.write (string+"\n")

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
	formatFilename = name (pdbFilename) + "-xyz-tmp-tinker.params"
	formatFile = open (formatFilename, "w")
	formatFile.write (pdbFilename + "\n")
	formatFile.write ("A" + "\n")
	formatFile.write (FORCE_FIELD + "\n")
	formatFile.close ()

	logFilename = createLogFilename (pdbFilename)

	cmm =  "pdbxyz < " + formatFilename  + " >> " + logFilename
	os.system (cmm)

	#cmm = "pdbxyz " + pdbFilename + " A " + FORCE_FIELD + " > " + logFilename
	#os.system ("pdbxyz " + pdbFilename + " A " + FORCE_FIELD + " > " + logFilename)
	clean ([formatFilename])
	xyzFilename = name (pdbFilename) + ".xyz"
	return xyzFilename

###############################################################################
# Create the temporal file including the parameters for 'analyze'program
###############################################################################
def createTinkerParametersFile (xyzFilename, FORCE_FIELD, typeOfAnalysis):
	proteinTinkerFilename = name (xyzFilename) + "anl--tmp-tinker.params"
	proteinFormatFile = open (proteinTinkerFilename, "w")

	proteinFormatFile.write (xyzFilename + "\n")
	proteinFormatFile.write (FORCE_FIELD + "\n")
	if typeOfAnalysis == "--potential":
		proteinFormatFile.write ("E")
	elif typeOfAnalysis == "--dipole":
		proteinFormatFile.write ("M")

	proteinFormatFile.close ()

	return proteinTinkerFilename

##################################################################
# Get the potential energy by preprocessing the output energy file
##################################################################
def getPotentialEnergy (logFilename):
	logFile = open (logFilename)
	for line in logFile:
		if "Total Potential Energy :" in line:
			value = line.split ()[4].strip() 
			value = value.replace ("D", "E")
			return value

	return "9999999999999999.999999999999"

##################################################################
# Get the dipole moment by preprocessing the output log file
##################################################################
def getDipoleMoment (logFilename):
	logFile = open (logFilename)
	for line in logFile:
		if " Dipole Moment Magnitude :" in line:
			value = line.split ()[4].strip() 
			return value

	return "9999999999999999.999999999999"

##################################################################
# RUN: Create xyz file, parameters file, and preprocess the energy output file
##################################################################
def run (pdbFilename, typeOfAnalysis):
	logFilename = createLogFilename (pdbFilename)
	xyzFilename = createXyzFromPdb (pdbFilename, FORCE_FIELD, logFilename)

	parametersFilename = createTinkerParametersFile (xyzFilename, FORCE_FIELD, typeOfAnalysis)

	os.system ("analyze < " + parametersFilename  + " >> " + logFilename )

	value = ""
	if typeOfAnalysis == "--potential":
		value = getPotentialEnergy (logFilename)
	elif typeOfAnalysis == "--dipole":
		value = getDipoleMoment (logFilename)

	clean ([logFilename, xyzFilename, parametersFilename, name(pdbFilename)+".seq"])
	return value

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	typeOfAnalysis = sys.argv[1]
	pdbFilename = sys.argv[2]

	print run (pdbFilename, typeOfAnalysis)


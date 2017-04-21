#!/usr/bin/python

USAGE  = "Calculate the number of native contacts using the MMSTB tool set\n"
USAGE += "USAGE: native-contacts-mm.py <PDB reference>  <PDB structure> [outputFilename]\n"

import os, sys 

###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> NatCont: "
	if type (message) != list:
		string +=  str (message)
	else:
		for i in message: string += str (i) + " "
	
	sys.stderr.write (string+"\n")

###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	import subprocess
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
	return value

############## CHECK ARGUMENTS ######################
if len (sys.argv) < 3:
	print USAGE
	sys.exit (0)
elif len (sys.argv) == 4:
	sys.stdout = open (sys.argv[3], "w")

############## MAIN #################################
pdbReference = sys.argv [1]
pdbFilename = sys.argv [2]

currentDir = os.getcwd()

log ([pdbReference, pdbFilename, currentDir])

strvalue = runProgram (["contact.pl", pdbReference, pdbFilename], currentDir)

value = strvalue.split ()[0]

print value


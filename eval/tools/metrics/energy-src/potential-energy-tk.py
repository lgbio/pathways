#!/usr/bin/python

USAGE = "Calculate the potential energy of a PDB using tinker functions"

import sys, os
###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	import subprocess
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
	return value

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

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv[1]

	outputDir = os.getcwd ()

	potentialEnergy = float ( runProgram (["energy", "--potential", pdbFilename], outputDir).strip())
	print potentialEnergy


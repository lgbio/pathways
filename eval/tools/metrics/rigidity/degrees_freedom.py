#!/usr/bin/python

"""
LOG:
	2016/12/16:	Added handling of standard error messages
"""

USAGE = "Calculates the Degrees of Freedom of a protein using FIRST rigity analysis program "

import sys, os
###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	import subprocess
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE, stderr=sys.stderr).communicate()[0]
	return value

###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> DegFree: "
	if type (message) != list:
		string +=  str (message)
	else:
		for i in message: string += str (i) + " "
	
	sys.stderr.write (string+"\n")

#-------------------------------------------------------------
# Define the  output for errors and log messages
#-------------------------------------------------------------
def defineMessagesOutput ():
	stderr = os.getenv ("EVAL_STDERR")
	if stderr == None:
		sys.stderr = sys.stdout
	else:
		sys.stderr = open (stderr, "a")

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":
	defineMessagesOutput ()

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv[1]

	outputDir = os.getcwd ()

	cmm = "rigidity-analysis.py %s %s" % ("--DegreesOfFreedom", pdbFilename)
	log (cmm)
	value = float ( runProgram (["rigidity-analysis.py", "--DegreesOfFreedom", pdbFilename], outputDir).strip())
	print value


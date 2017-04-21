#!/usr/bin/python

"""
LOG:
	2016/12/16:	Added handling of standard error messages
"""

USAGE = "Calculate the proportion of amino acids in correct positions"

import sys, os
###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	import subprocess
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE, stderr=sys.stderr).communicate()[0]
	return value

#-------------------------------------------------------------
# Print a message to the error output stream
#-------------------------------------------------------------
def log (message):
	string=">>> SecStrCrr: "
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

	pdbTarget = sys.argv[1]
	pdbReference = sys.argv[2]

	outputDir = os.getcwd ()

	value = float ( runProgram (["secondary-structures", "-correct", pdbTarget, pdbReference], outputDir).strip())
	print value


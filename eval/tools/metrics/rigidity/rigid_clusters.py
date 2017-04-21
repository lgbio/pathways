#!/usr/bin/python

USAGE = "Calculates the RIGID CLUSTERSt of a protein using FIRST rigity analysis program "

import sys, os
###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	import subprocess
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
	return value

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv[1]

	outputDir = os.getcwd ()

	value = float ( runProgram (["rigidity-analysis.py", "--RigidClusters", pdbFilename], outputDir).strip())
	print value


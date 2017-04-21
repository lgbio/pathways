#!/usr/bin/python

USAGE  = "Calculate the number of hydrogen bonds using the hbpluss tool\n"
USAGE += "USAGE: hbonds-lg.py <PDB structure> [outputFilename]\n"

import os, sys

############## MAIN #################################
def main (args):
	if len (sys.argv) < 2:
		print USAGE
		sys.exit (0)
	elif len (sys.argv) == 3:
		sys.stdout = open (sys.argv[2], "w")

	pdbFilename     = args [1]
	pdbName         = os.path.basename (pdbFilename).split(".")[0]
	logFilename     = pdbName + ".log"
	hbondsFilename  = pdbName + ".hb2"

	command = "hbplus " + pdbFilename + " >& " + logFilename 
	log (command)
	os.system (command)

	resultsOfHbonds = open (logFilename).readlines()[-1].split ()[0]

	os.remove (logFilename)
	os.remove (hbondsFilename)

	print resultsOfHbonds

#-------------------------------------------------------------
# Print a message to the error output stream
#-------------------------------------------------------------
def log (message):
	string=">>> HBonds: "
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

#-------------------------------------------------------------
# MAIN
#-------------------------------------------------------------
if __name__ == "__main__":
	defineMessagesOutput ()
	main (sys.argv)

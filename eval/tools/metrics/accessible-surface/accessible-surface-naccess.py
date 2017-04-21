#!/usr/bin/python

"""
AUTHOR: Luis Garreta
DATE: May 07/2010
	  modified Dic 14/2017

"""
USAGE  = "Return the total  polar and non-polar superfice accessible to solvent using 'naccess' program\n"
USAGE += "USAGE:  sas-naccess.py <pdb filename>\n"
USAGE += "INPUT:  the pdb filename\n"
USAGE += "OUTPUT: Two value: the all-pollar SAS and the non-polar SAS"

import sys, os

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

###############################################################################
# 
###############################################################################
def run (pdbFilename):
	os.environ ["NACCESS_HOME"] = os.getenv ("EVAL_HOME") + "/tools/programs/naccess"
	pathname, filename = os.path.split (pdbFilename)
	stemName = name (pdbFilename)

	# Temporal files created by the naccess program
	rsaFilename = pathname + "/" + filename.split (".")[0] + ".rsa"
	logFilename = pathname + "/" + filename.split (".")[0] + ".log"
	asaFilename = pathname + "/" + filename.split (".")[0] + ".asa"
	inputFilename = pathname + "/" + filename.split (".")[0] + ".input"

	if pathname == "": pathname = "."

	cmm =  "naccess " + pdbFilename + " >> /dev/null " 
	log (cmm)
	#log ([logFilename, rsaFilename, asaFilename, inputFilename])
	os.system (cmm)

	# Process rsa file to get the total ASA
	rsaFile = open (rsaFilename)
	listOfLines = rsaFile.readlines()
	totalValues = listOfLines [-1].split ()
	totalSAS = totalValues [1]

	# Remove temporal files
	clean ([logFilename, rsaFilename, asaFilename, inputFilename])

	return totalSAS

#-------------------------------------------------------------
# Print a message to the error output stream
#-------------------------------------------------------------
def log (message):
	string=">>> ASA: "
	if (type (message) != list): 
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

	sas = run (pdbFilename)
	print sas 


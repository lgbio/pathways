#!/usr/bin/python

"""
Return the total void volume in a protein. It uses the 'avp' program.
It uses a default probe of 0.5 Angstrom

AUTHOR: Luis Garreta
DATE: Dic 06/2010

"""
USAGE  = "Return the total void volume in a protein. It uses the 'avp' program\n"
USAGE += "USAGE	:  voids-avp.py <pdb filename>\n"
USAGE += "INPUT	:  the pdb filename\n"
USAGE += "OUTPUT: The total void volume of the protein"

import sys, os
###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=">>> Voids: "
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

###############################################################################
# 
###############################################################################
def run (pdbFilename):
	stemName = name (pdbFilename)

	logFilename = stemName + ".log"
	outFilename = stemName + ".out"
	
	params = {"IN": pdbFilename, "OUT": outFilename, "LOG":logFilename}
	command = "avp -r %(IN)s %(OUT)s &> %(LOG)s" % params
	os.system (command)

	outFile = open (outFilename)
	listOfLines = outFile.readlines()

	totalValues = listOfLines [-2].split ()

	totalVoids = totalValues [-1]

	clean ([outFilename, logFilename])

	return totalVoids

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv[1]

	out = run (pdbFilename)
	print out 


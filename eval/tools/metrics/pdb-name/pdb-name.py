#!/usr/bin/python

"""
Return the name of a PDB file

AUTHOR: Luis Garreta
DATE: March 12/2012
"""
USAGE  = "Return the name of a PDB file\n"
USAGE += "USAGE	:  pdb-name.py <pdb filename>\n"
USAGE += "INPUT	:  the pdb filename\n"
USAGE += "OUTPUT: The name of the pdb file "

import sys, os
###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
	baseName = os.path.basename (namefile)
	newNamefile = baseName.split(".") [0]

	return newNamefile

###############################################################################
# 
###############################################################################
def run (pdbFilename):
	return name (pdbFilename)

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


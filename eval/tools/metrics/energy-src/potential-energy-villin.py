#!/usr/bin/python

USAGE = "Calculate the potential energy of a PDB using tinker functions"

import sys, os
ENERGIES_FILE = "/home/data/phd/pathways/villin/PROJ3036/RUN4/all/all-energies-run4-clone0.values"
###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
from numpy import *
def getEnergy (pdbFilename):
	## get time
	pdbFile = open (pdbFilename)
	pdbFile.next()
	timeLine =  pdbFile.next()
	timeStr = timeLine.split ("t=")[1].split (".")[0]
	time = (int (float (timeStr))) / 50

	## get energy
	energiesArray = loadtxt (ENERGIES_FILE, dtype={'names':('time', 'energy'), 'formats':('i4','f4')})

	energy = energiesArray [time][1]

	return energy

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	pdbFilename = sys.argv[1]

	energy = getEnergy (pdbFilename)

	print energy




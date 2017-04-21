#!/usr/bin/python

import sys

##################################################################
# Get energy terms from a Gromas energy file (.edr)
##################################################################

def getGromacsEnergyTerms (energyFilename):
	energyFile = open (energyFilename)
	energyLines = energyFile.readlines()
	listOfValues = energyLines [-1].strip().split()

	potentialEnergy = listOfValues [1]
	kineticEnergy = listOfValues [2]
	totalEnergy = listOfValues [3]
	temperature = listOfValues [4]
	pressure = listOfValues [5]

	return potentialEnergy, kineticEnergy, totalEnergy, temperature, pressure

##################################################################
# MAIN
##################################################################
if __name__ == "__main__":

	if len (sys.argv) < 2: 
		print USAGE
		sys.exit (0)

	energyFilename = sys.argv[1].split (".")[0] + ".xvg"

	print getGromacsEnergyTerms (energyFilename)


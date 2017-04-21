#!/usr/bin/python

import os, sys

USAGE  = "Returns the potential energy of a PDB conformation using gromacs g_energy"
USAGE += "potential-energy.py <pdb file>"

def getForceField ():
	forceField = os.getenv ("PHD_FORCE_FIELD")
	if forceField == "AMBER96":
		return "amber96"
	elif forceField == "OPLSAA":
		return "oplsaa"
	elif forceField == "GROMOS":
		return "G43a1"
###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
	string=""
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

###############################################################################
# Create a gromacs structure (.tpr) file from a protein file (.pdb)
###############################################################################
def potentialEnergy (pdbName):
	tmpName = name (pdbName) + "-tmp-energy"
	forceField = getForceField ()

	command = "pdb2gmx -f %s -o %s -p %s -i %s -ff %s -ignh -missing >> %s.log 2>&1" % (pdbName, tmpName, tmpName, tmpName, forceField, tmpName)
	log ([command])
	os.system (command)

	command = "editconf -f %s -o %s -d 0.5 >> %s.log 2>&1" % (tmpName, tmpName, tmpName)
	log ([command])
	os.system (command)

	command = "genbox -cp %s -cs -p %s -o %s >> %s.log 2>&1" % (tmpName, tmpName, tmpName, tmpName)
	log ([command])
	os.system (command)

	command="ln -s $PHD_SCRIPTS_MEASURES_PATH/em-energy.mdp %s.mdp >> %s.log 2>&1" % (tmpName, tmpName)
	log ([command])
	os.system (command)

	command = "grompp -f %s -c %s -p %s -po %s -o %s >> %s.log 2>&1 " % (tmpName, tmpName, tmpName, tmpName, tmpName, tmpName) 
	log ([command])
	os.system (command)

	command = "mdrun -s %s -o %s -c %s -e %s -cpo %s -g %s >> %s.log 2>&1" % (tmpName, tmpName, tmpName, tmpName, tmpName, tmpName, tmpName)
	log ([command])
	os.system (command)

	command = "echo Potential | g_energy -f %s -o %s >> %s.log 2>&1" % (tmpName, tmpName, tmpName)
	log ([command])
	os.system (command)

	lastLine = open (tmpName + ".xvg").readlines()[-1]
	potentialEnergyValue = lastLine.split()[-1]

	command = "rm *%s*" % os.path.basename (tmpName)  
	log ([command])
	os.system (command)

	return potentialEnergyValue

###############################################################################
###############################################################################

if __name__ == "__main__":

	if len (sys.argv) < 2 :
		##print USAGE
		sys.exit (0)

	pdb = sys.argv [1]

	print potentialEnergy (pdb)

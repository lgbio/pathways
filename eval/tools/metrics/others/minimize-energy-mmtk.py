#!/usr/bin/python 

"""
MMTK has two energy minimizers using different algorithms: steepest descent
and conjugate gradient. The former has the advantage of always moving towards 
the minimum that is closest to the starting point and is therefore ideal for 
removing bad contacts in a unreasonably high energy configuration.
"""

from MMTK import *
from MMTK.Proteins import *
from MMTK.ForceFields import Amber94ForceField
from MMTK.Minimization import SteepestDescentMinimizer

#####################################################################
#####################################################################
def minimizeEnergy (inputProteinName, outputProteinName):
	universe = InfiniteUniverse(Amber94ForceField())
	universe.protein = Protein (inputProteinName)
	minimizer = SteepestDescentMinimizer(universe)
	minimizer(steps = 10)

	universe.protein.writeToFile (outputProteinName)

#####################################################################
import sys
if __name__ == "__main__":
	if len (sys.argv) < 3:
		print USAGE
		sys.exit (0)


	inputProtein = sys.argv [1]
	outputProtein = sys.argv [2]

	minimizeEnergy (inputProtein, outputProtein)


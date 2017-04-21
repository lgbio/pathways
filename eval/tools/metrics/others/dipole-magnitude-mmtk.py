#!/usr/bin/python

## Return the vector magnitude of the dipole moment

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
import sys, math

pdb = sys.argv [1]
test = Protein(pdb)
universe = InfiniteUniverse(Amber94ForceField())
universe.addObject(test)
vector = universe.dipole()
magnitude = math.sqrt (pow (vector [0],2) + pow (vector [1],2) + pow (vector [2],2))

print magnitude

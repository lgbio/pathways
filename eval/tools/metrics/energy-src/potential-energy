#!/usr/bin/python
from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
import sys

pdb = sys.argv [1]
test = Protein(pdb)
universe = InfiniteUniverse(Amber99ForceField())
universe.addObject(test)
print universe.energy()


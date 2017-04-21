from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
import sys

pdb = sys.argv [1]
test = Protein(pdb)
universe = InfiniteUniverse(Amber94ForceField())
universe.addObject(test)

print universe.energy()
for i in universe.energyTerms ():
	print i, ":", 
	print universe.energyTerms()[i]

#print "angular moment:", test.angularMomentum()
print "charge :", universe.charge()
print "degrees of freedom :", universe.degreesOfFreedom()
print "dipole moment:", universe.dipole()
#print "kinetic energy :", test.kineticEnergy()
print "mas :", universe.mass()
#print "momentum :", g.momentum()
print " GroupOfAtoms numberOfAtoms:", universe.numberOfAtoms ()
#print " temperature:", universe.temperature ()

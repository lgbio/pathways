#!/usr/bin/python
"""
The Python script below will calculate the radius of gyration for the
assembly of all molecules specified in a PDB file (typically the
		asymmetric unit). To run it, you need the Molecular Modelling
Toolkit, available from http://dirac.cnrs-orleans.fr/MMTK/
"""

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.PDBMoleculeFactory import PDBMoleculeFactory
from Scientific import N
import sys

conf = PDBConfiguration(sys.argv[1])
factory = PDBMoleculeFactory(conf)
molecules = Collection(factory.retrieveMolecules())

def radiusOfGyration(m):
	natoms = m.numberOfAtoms()
	center = sum((atom.position() for atom in m.atomList()),
			Vector(0., 0., 0.)) / natoms
	sum_r = sum(((atom.position()-center).length()**2
		for atom in m.atomList()))
	return N.sqrt(sum_r/natoms)

print radiusOfGyration(molecules)/Units.Ang

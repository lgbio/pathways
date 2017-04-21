#!/usr/bin/python

USAGE  = "Calculate the number of native contacts using the biopython lib\n"
USAGE += "USAGE: native-contacts <PDB reference>  <PDB structure>\n"
"""
NOTES: The calculation only take into account 3 far from residues and a 
       minimal distance of 7A
"""

import Bio.PDB
import numpy
import sys

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two, min, max) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.int)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            distance = calc_residue_dist(residue_one, residue_two)
            if abs (row-col) >= min:
                if distance <= max:
				    answer[row, col] = 1
                else:
				    answer[row, col] = 0

    return answer

def calc_nativeContacts (matRef, matCmp):
	mat = matRef & matCmp
	c = 0
	i = 0
	n = len (mat)
	for i in range (0,n):
		for j in range (i+1,n):
				if mat [i][j] == 1:
					c+=1
	return c
#####################################################################
## Main	
#####################################################################

if len (sys.argv) < 3:
	print USAGE
	sys.exit (0)

pdb_reference = sys.argv [1]
pdb_comparison = sys.argv [2]

structureRef = Bio.PDB.PDBParser().get_structure ("Ref", open (pdb_reference))
structureCmp = Bio.PDB.PDBParser().get_structure ("Cmp", open (pdb_comparison))

modelRef = structureRef[0]
modelCmp = structureCmp[0]

matrixRef = calc_dist_matrix (modelRef.child_list[0], modelRef.child_list[0], 3, 7)
matrixCmp = calc_dist_matrix (modelCmp.child_list[0], modelCmp.child_list[0], 3, 7)

nativeContacts = calc_nativeContacts (matrixRef, matrixCmp)

print nativeContacts

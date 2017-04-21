#!/usr/bin/Rscript

# Calculate the RMSD of two proteins.
# It fits the second one --coordinte superposition.
# Proteins are entered as command line arguments
# USAGE: rmsd.R <proteinA> <proteinB> [outputFilename]

library (bio3d)

args <- commandArgs (TRUE)

proteinA <- args[1]; 
proteinB <- args[2]

#--- Obtain the stem filename of the compared protein (proteinB)
proteinBFilename = unlist (strsplit (proteinB, "\\."))[1]
rmsdFilename = paste (proteinBFilename, "", ".rmsd", sep="")
if (length (args) == 3)
	rmsdFilename <- args [3]

A <- read.pdb (proteinA, verbose=FALSE)
A.ind <- atom.select (A, elety="CA", verbose=FALSE)
A.ind

B <- read.pdb (proteinB, verbose=FALSE)
B.ind <- atom.select (B, elety="CA", verbose=FALSE)
B.ind


#--- Calculte the RMSD fitting the two proteins (coordinate superpostion)
r1 = rmsd (A$xyz[A.ind$xyz], B$xyz[B.ind$xyz], fit=TRUE)
write (r1, rmsdFilename)
cat (r1)

##--- Fit two PDBs
#xyz <- fit.xyz(fixed=A$xyz, mobile=B$xyz, fixed.inds=A.ind$xyz, mobile.inds=B.ind$xyz)

#-- Write out moved PDB
#movedNamefile = paste (proteinBFilename, "", "-moved.pdb", sep="")
#C <- B; C$xyz = xyz
#write.pdb(pdb=C, file = movedNamefile)


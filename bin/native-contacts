#!/usr/bin/Rscript

# Calculated the number of native contacts of a pdb structure
# USAGE: contacts <PDB reference file> <PDB file>

# LOG:
#    Dic042012: Change checking number of parameters

suppressPackageStartupMessages (library (bio3d))

args <- commandArgs (TRUE)

if (length (args) < 2) {
	print ("USAGE: native-contacts <Target> <Reference>")
	quit()
}

proteinA <- args[1] 
proteinB <- args[2]

# Protein A
sink ("/dev/null")
pdbA <- read.pdb (proteinA, verbose=FALSE)
indA <- atom.select (pdbA, "calpha", verbose=FALSE)

# Protein B
pdbB <- read.pdb (proteinB, verbose=FALSE)
indB <- atom.select (pdbB, "calpha", verbose=FALSE)
sink(NULL)
contA <- cmap (pdbA$xyz [indA$xyz], dcut=7, scut=3)

contB <- cmap (pdbB$xyz [indB$xyz], dcut=7, scut=3)

prod <- contA * contB
sum <- sum (prod, na.rm=TRUE)

write (sum, file="")

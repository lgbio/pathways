#!/usr/bin/Rscript

library (bio3d)
library (parallel)
library (cluster)

source ("libs/createDir.R")
source ("libs/splitFilesToBins.R")

nCPUS = 4
USAGE = "USAGE: medoids.R <input dir> <output dir> <size chunks>\n"
options (width=400)
tmpDir = "/dev/shm"
#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function (args) {
	#args = c("shm/1FCA1", "1000")
	if (length (args) < 1) {
		cat (USAGE)
		quit ()
	}
	inputPathname  = args [1]
	sizeMedoids    = strtoi (args [2])
	if (length (args) > 2) 
		outputDir      = args [3]
	else
		outputDir      = getwd ()

	cat (">>> Splitting Files into bins...\n")
	binDirList = splitFilesToBins (inputPathname, sizeMedoids, tmpDir)

	medoidsDir =  sprintf ("%s/%s", outputDir, "medoids")
	createDir (medoidsDir)
	cat (">>> Creating medoids...\n")
	for (binDir in binDirList) {
		medoid = medoidsFromBin (binDir)
		system (sprintf ("cp %s %s/%s", medoid, outputDir, "medoids"))
	}
}

#--------------------------------------------------------------
#--------------------------------------------------------------
medoidsFromBin <- function (binDir) {
	outputGroups   = paste (binDir, ".groups", sep="")
	outputMedoids  = paste (binDir, ".medoids", sep="")

	cat (">>> Loading pdb filenames...\n")
	pdbObjects <<- getPDBFiles (binDir)

	cat ("\n>>> Clustering PDBs...\n")
	results = clusteringPDBs (pdbObjects)

	cat ("\n>>> Writing Results...\n")
	write.table (file=outputGroups, results$groups)
	write.table (file=outputMedoids, results$medoids)

	medoidFilename = list.files (binDir) [results$medoids]
	return (sprintf ("%s/%s", binDir, medoidFilename))
}

#--------------------------------------------------------------
clusteringPDBs <- function (pdbObjects) {
	cat ("   >>> Calculating RMSD Distances...\n")
	rmsdDistances <<- pairwiseDistancesRMSDs (pdbObjects)

	pamPDBs <<- pam (rmsdDistances, 1, diss=T)
	groups  <- pamPDBs$clustering
	medoids <- pamPDBs$medoids

	return (list(groups=groups,medoids=medoids))
}


#--------------------------------------------------------------
# Calculate pairwise RMSDs
#--------------------------------------------------------------
pairwiseDistancesRMSDs <- function (pdbObjects) {
	pdbs = pdbObjects$pdbs
	firstPdb = pdbs [[1]]
	n = length (pdbObjects$pdbs)

	# Calculate matrix of coordinates xyz
	xyzMatrix <- matrix (firstPdb$xyz,nrow=1)
	for (pdb in pdbs [2:n]) {
		xyzMatrix = rbind (xyzMatrix, n1=pdb$xyz)
	}
	rownames (xyzMatrix) <- names (pdbObjects$pdbs)

	# Calculate RMSDs
	# Trajectory Frame Superposition on Calpha atoms
	ca.inds <- atom.select(firstPdb, "calpha")
	xyz <- fit.xyz (fixed = firstPdb$xyz, mobile = xyzMatrix, 
		fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz, ncore=nCPUS)

	#rmsdDistances <- as.dist (rmsf (xyz, average=T))
	rmsdDistances <- as.dist (rmsd (xyz, a.inds=ca.inds$xyz, ncore=nCPUS))
	#rd <- rmsd (xyz[1,], xyz, a.inds=ca.inds$xyz, ncore=nCPUS)
	#rownames (rd) = rownames (xyzMatrix)

	return (rmsdDistances)
}
#--------------------------------------------------------------
# Calculate matrix of XYZ
#--------------------------------------------------------------
calcXYZMatrix <- function (pdbObjects, native) {
	n = length (pdbObjects)

	xyzMatrix = matrix (native$xyz,nrow=1)
	for (pdb in pdbObjects [2:n]) {
		xyzMatrix = rbind (xyzMatrix, n1=pdb$xyz)
	}
	rownames (xyzMatrix) <- names (pdbObjects)
	return (xyzMatrix)
}

#--------------------------------------------------------------
# Get RMSDs between each pair of protein structures
#--------------------------------------------------------------
calcRMSDsPairwise <- function (native, xyzMatrix) {
	# Trajectory Frame Superposition on Calpha atoms
	ca.inds <<- atom.select(native, "calpha")
	xyz <<- fit.xyz (fixed = native$xyz, mobile = xyzMatrix, 
		fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz, ncore=1)

	rd <- rmsd (xyz, a.inds=ca.inds$xyz, ncore=nCPUS)
	rownames (rd) = rownames (xyzMatrix)

	return (rd)
}

#--------------------------------------------------------------
# Load pdb files to pdb objects
#--------------------------------------------------------------
getPDBFiles <- function (pathname) {
	cat ("\n\n>>>path ", pathname, "\n")
	# Extracts to an pathname 
	if (grepl ("gz", pathname)==T) {
		stemName = strsplit (pathname, split="[.]")[[1]][1]
		untar (pathname, compressed=T, exdir=stemName)
		inputDir = stemName 
	} else 
		inputDir = pathname

	pdbNames     = list.files (inputDir)
	pdbNamesFull = sapply (pdbNames, function (x) paste (inputDir, x, sep="/"))
	n = length (pdbNamesFull)
	native = pdbNamesFull [[n]]

	# Load PDB Objects
	nativeObject <<- read.pdb2 (native)
	pdbObjects <<- mclapply (X=pdbNamesFull, FUN=read.pdb2, mc.cores=nCPUS )
	
	return (list (n=n, native=nativeObject, pdbs=pdbObjects))
}

#--------------------------------------------------------------
# Print to text file the input object
#--------------------------------------------------------------
log <- function (object, dwfilename) {
	sink (filename)
	print (object)
	sink()
}

#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------

args = commandArgs (TRUE)
main (args)

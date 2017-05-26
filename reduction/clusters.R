#!/usr/bin/Rscript

library (bio3d)
library (parallel)

#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function () {
	args = commandArgs (TRUE)
	#inputDir = args [1]
	inputDir = "in"
	inputFullFiles = paste (inputDir, list.files (inputDir), sep="/")

	print (">>>")
	print ( inputFullFiles)
	mcmapply (clusteringProcess, inputFullFiles, mc.cores=4)
}

#--------------------------------------------------------------
clusteringProcess <- function (pathname) {
	#pathname  = args [1]  # Must be as in/1FA01.tgz
	outFile =  gsub (".tgz", ".tbl", gsub ("in/", "out/", pathname)) 

	cat ("Getting pdb filenames...\n")
	files <<- getPDBFiles (pathname)

	cat ("Calculating XYZ Matrix...\n")
	xyzMatrix <<- calcXYZMatrix (files$pdbs)

	cat ("Calculating RMSDs..\n")
	rmsdDistances <<- calcRMSDsPairwise (files$native, xyzMatrix)

	cat ("Running clustering...")
	groups <<- runClustering (rmsdDistances)

	write.table (file=outFile, groups)
}

#--------------------------------------------------------------
# Calculate matrix of XYZ
#--------------------------------------------------------------
calcXYZMatrix <- function (pdbNames) {
	n = length (pdbNames)

	firstPdb <<- read.pdb2 (pdbNames[1])
	#ca.inds <- atom.select(firstPdb, "calpha")
	#firstPdb <- trim (firstPdb, ca.inds)
	xyzMatrix = matrix (firstPdb$xyz,nrow=1)

	pdbs <<- mclapply (pdbNames, read.pdb2)
	#print (pdbs)
	for (pdb in pdbs [2:n]) {
		#pdb = read.pdb2 (pdbName)
		#ca.inds <- atom.select(pdb, "calpha")
		#pdb <- trim (pdb, ca.inds)
		#xyzMatrix = rbind (xyzMatrix, n1=pdb$xyz)
		xyzMatrix = rbind (xyzMatrix, n1=pdb$xyz)
	}
	rownames (xyzMatrix) <- names (pdbNames)
	return (xyzMatrix)
}

#--------------------------------------------------------------
# Get RMSDs between each pair of protein structures
#--------------------------------------------------------------
runClustering <- function (rmsdDistances) {
	ds <- dist (rmsdDistances)
	hs <- hclust (ds)
	groups <- cutree (hs, k=4)

	return (groups)
}

#--------------------------------------------------------------
# Get RMSDs between each pair of protein structures
#--------------------------------------------------------------
calcRMSDsPairwise <- function (native, xyzMatrix) {
	# Trajectory Frame Superposition on Calpha atoms
	native <- read.pdb2 (native)
	ca.inds <<- atom.select(native, "calpha")
	#native <<- trim (native, ca.inds)
	xyz <<- fit.xyz (fixed = native$xyz, mobile = xyzMatrix, 
		fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz, ncore=1)
	# Root Mean Square Deviation (RMSD)

	#rownames (xyz) = rownames (xyzMatrix)

	rd <- rmsd (xyz, ncore=1)
	rownames (rd) = rownames (xyzMatrix)

	#plot(rd, typ = "l", ylab = "RMSD", xlab = "Frame No.")
		#points (lowess(rd), typ = "l", col = "red", lty = 2, lwd = 2)

	return (rd)
}

#--------------------------------------------------------------
# Get the PDB files from either a compressed file or a dir
#--------------------------------------------------------------
getPDBFiles <- function (pathname) {
	# Extracts to an inputDir 
	if (grepl ("gz", pathname)==T) {
		stemName = strsplit (pathname, split="[.]")[[1]][1]
		inputDir = stemName 
		untar (pathname, compressed=T, exdir=inputDir)
	} else 
		inputDir = pathname

	inputFiles = list.files (inputDir)
	inputFilesFull = sapply (inputFiles, function (x) paste (inputDir, x, sep="/"))
	#inputFilesFull = mclapply (inputFiles, function (x) paste (inputDir, x, sep="/"))
	n = length (inputFilesFull)
	native = inputFilesFull [[n]]
	outfile = paste (inputDir, ".pdf", sep="")

	return (list (n=n, native=native, pdbs=inputFilesFull, outfile=outfile))
}

#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------
main ()

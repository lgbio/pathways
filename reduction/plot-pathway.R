#!/usr/bin/Rscript

library (bio3d)
library (parallel)

MC_CORES=4

# It takes a pathway and reduces by three methods:
# static, random, mediodes, clustering.
# Plots are created for each result.

#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function () {
	args = commandArgs (TRUE)
	pathname  = args [1]

	# Extract or load filename to calculate RMSDs
	files = getPDBFiles (pathname)
	
	rmsdValues = parallerRmsdPathway (files$n, files$native, files$pdbs)

	# Plot for the original pathway values
	plotPathway (rmsdValues, files$outfile)
	#plot (rmsdValues, type="l", lty="solid")
}

#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
parallerRmsdPathway <- function (n, native, pdbs) {
	rmsdValues = mclapply (pdbs, calculateRMSD, native, mc.preschedule=F)
	return (rmsdValues)
}
	
rmsdPathway <- function (n, native, pdbs) {
	rmsdValues = c()
	for (i in 1:n) {
		rmsd = calculateRMSD (pdbs [i], native)
		rmsdValues = append (rmsdValues, rmsd)
	}
	return (rmsdValues)
}

#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
calculateRMSD <- function (pdbNameTarget, pdbNameRef) {
	print (sprintf ("%s, %s\n", pdbNameTarget, pdbNameRef))
	#--- Obtain the stem pathname of the compared protein (pdbNameRef)
	proteinTargetFilename = unlist (strsplit (pdbNameTarget, "\\."))[1]

	target <- read.pdb2 (pdbNameTarget, rm.alt=FALSE, verbose=FALSE)
	reference <- read.pdb2 (pdbNameRef, rm.alt=FALSE, verbose=FALSE)
	lenAtomsTarget <- length (target$atom[,1])
	lenAtomsReference <- length (reference$atom[,1])

	targetCAs <- atom.select (target, elety="CA", verbose=FALSE)
	referenceCAs <- atom.select (reference, elety="CA", verbose=FALSE)

	#--- Calculte the RMSD fitting the two proteins (coordinate superpostion)
	r1 = rmsd (target$xyz[targetCAs$xyz], reference$xyz[referenceCAs$xyz], fit=TRUE)

	return  (r1)
}


#-------------------------------------------------------------
# Creata a XY plot from the RMSD values of each conformation
#-------------------------------------------------------------
plotPathway <- function (rmsdValues, outfile) {
	pdf (outfile, width=2.8, height=2)
		# Outer margins, scale fonts, pos labels (
		par (mar=c(3,3,1,1), cex=0.7, mgp = c(2,1,0)) 
		n = length(rmsdValues)
		time = 1:(n-1)
	 	plot (time, rmsdValues[1:(n-1)], xlab="Steps", ylab="RMSD", type="l", lty="solid", col=1)
		#legend (x="top", legend=c("Original","Averaged"), lty=c(2,1), col=c(1,2))

		par (new=T)
	 	#lines (avrPropValues, col="red", lty=1)
	dev.off()
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
	n = length (inputFilesFull)
	native = inputFilesFull [[n]]
	outfile = paste (inputDir, ".pdf", sep="")

	return (list (n=n, native=native, pdbs=inputFilesFull, outfile=outfile))
}

#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------
main ()

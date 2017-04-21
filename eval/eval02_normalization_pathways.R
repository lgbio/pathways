#!/usr/bin/Rscript

# Normalize multiple pathways
# Normalize a set of values from an input file by unitization
# na = (a - min(a)) / (max(a)-min(a))


options (width=400)
args = commandArgs (TRUE)

main <- function () {
	dirName <- args[1]
	normalizeMultiplePathways (dirName)
}

normalizeMultiplePathways <- function (dirName) {
	inputFiles = list.files (dirName)	
	fullInputFiles = sapply (inputFiles, function (x) paste (dirName, x, sep="/"))
	sapply (fullInputFiles, normalizeSinglePathway  )
}	

normalizeSinglePathway <- function (dataFilename){
	print (paste (">>>", dataFilename))
	#dataFilename <<- "sample-10-prop.values"
	baseFilename <<- strsplit (dataFilename, split="\\.")[[1]][1]
	newFilename  <<- sprintf ("%s-norm.values", baseFilename)

	values <<- read.table (dataFilename, header=T, row.names=1)

	normalizedValues = apply (values, 2, normalization)
	formatedValues   = apply (normalizedValues, 2, formatValue)
	write.table (file=newFilename, formatedValues, sep=" ", quote=F)
}

formatValue <- function (value)  {
	return (format (value, digits=6))
}


#--------------------------------------------------------------------
#--------------------------------------------------------------------
normalization <- function (x){
	minX = min (x)
	maxX = max (x)

	if (maxX == minX) return (0)

	newX = (x - minX) / (maxX - minX)
	return (newX)
}


#--------------------------------------------------------------------
#--------------------------------------------------------------------
main()

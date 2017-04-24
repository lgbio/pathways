#!/usr/bin/Rscript

### GOAL: 
###    Make a clustering from getPrincipalComponents scores and classify a pathway into these
###    groups using the loadinga from getPrincipalComponentss and the evaluated properties.
###    Create two graphics: Components&Groups and Pathway classification
### LOG: 

	###    Aug   19/15: Modified Plot Clustering to Write Y labels
	###    Aug   16/15: Clustering Run with the full dataset but data changes (not analyzed)

	###    Oct   25/14: Fixed with the framework function
	###    May   13/14: Added boxplot characteristics on levels (boxplotsComponentsOnGroups)
	###    May   08/14: Added curve line in boxplots for groups
	###    May   08/14: Fixed getRotatedPrincipalComponents (old version, no loop)
	###    Nov   10/13: Added own PCA function (Preliminary calculated with this PC)
	###    Nov   09/13: Main as a first function

loadConstants <- function () {
	AVR              <<- 5           # Average of 5 (2 left, 1 center, 2 right) smoothing noisy values
	CLUSTER_DISTANCE <<- "euclidean" # euclidean distance for clustering
	CLUSTER_LINKAGE  <<- "complete"  # Complete linkage method for hierarchical clustering
	N_COMPONENTS     <<- 3
	N_CLUSTERS       <<- 4
	BAD_PROPERTIES   <<- c("PE", "SR", "DM", "RC", "VD")  # Removed from the full values table
	datasetName      <<- sprintf ("inputs/parasol-dataset-avr%s-selected.values", AVR) # Name of the file dataset
	LABEL            <<- name ("-avr")   # Labels for output files average 5, variance, rotated, cross val 20%t
	outputDir        <<- name ("out", LABEL)
	DILUTION_FACTOR  <<- 0.8 	   ## dilutionFactor    # Porcentaje of training
}

options (width=300)
##################################################################
# MAIN 
# Run measure analysis by selecting various datasets for training and test (cross validattion)
##################################################################
main <- function () {
	loadConstants ()

	datasetName      = sprintf ("inputs/parasol-dataset-avr%s-selected.values", AVR) # Name of the file dataset
	outputDir       = name ("output-","pca-analysis")
	clusterFilename = name ("pca-cluster-parasol.RData")
	createOutputDir (outputDir)
	allValues       = read.table (datasetName)

	## Remove "noisy" measures, four distant plus PE 
	valuesCleaned     = cleanValues (allValues, BAD_PROPERTIES)
	pathwayNames      = getPathwayNames (valuesCleaned)

	## Label for output files average 5, variance, rotated, cross val 20%
	LABEL             = name ("-avr5-val", DILUTION_FACTOR)   
	dataSplitted      = splitDataset (pathwayNames, DILUTION_FACTOR, 1) # Return only one split
	valuesTraining    = getValuesPathways (valuesCleaned, dataSplitted$trainingNames)
	valuesTesting     = getValuesPathways (valuesCleaned, dataSplitted$testingNames)

	print (outputDir)
	setwd (outputDir)

	measureAnalysisSingleVersion (allValues, valuesTraining, valuesTesting, clusterFilename, LABEL)

	namesAxis   = c("Conformational Features", "Scores")
	namesPCs    = c("Stability","Compactness","Native-likeness")
	#namesGroups = c("Group 1", "Group 2", "Group 3", "Group 4")
	namesGroups = c("Group g1", "Group g4", "Group g2", "Group g3")  # For paper natural cluster dendr.
	labelFile   = "comp-"
	boxplotsComponentsOnGroups (clusterH, namesPCs, 
	                            namesGroups, namesAxis, labelFile)

	namesAxis   = c("Conformational Features", "Values")
	namesGroups = c("Unfolded", "Early Intermediate", "Late Intermediate", "Folded")
	labelFile   = "char-"
	boxplotsComponentsOnGroups (clusterH, namesPCs, 
	                            namesGroups, namesAxis, labelFile)

	setwd ("..")

	#measureAnalysisResamplingVersion (20)
}

###################################################################################
# Single evaluation
# Make clustering with all the values from the dataset
###################################################################################
measureAnalysisSingleVersion <- function (allValues, valuesTraining, valuesTesting, clusterFilename, label) {
	## Preliminary MDS to visualize an underlying struture

	BAD_PROPERTIES   <<- c("PE", "SR", "DM", "RC", "VD")  # Removed from the full values table
	MDSAnalysis (allValues, BAD_PROPERTIES, "-pre")
	quit()

	## Preliminary PCA to filter data. First, full measures, second without the four distant
	#pcaPre = preliminaryPCA (allValues)

	## Fina boxplot(F~classes, scc, names=c("Unfolded", "Early Interm.", "Late Interm.","Folded"), col=c(8,2,3,4),xlab="Folding Levels", ylab="Score", main="ICF Score")ll PCA to filter data
	analysis  = measureAnalysis (valuesTraining, valuesTesting, clusterFilename, 3, label)

	clusterH   <<- analysis$cluster 
	plotClusteringInfo (clusterH, N_CLUSTERS, "" )
	boxplotsMatrix (clusterH$scores, N_CLUSTERS, individualFiles=F, "components")
	boxplotsMatrix (clusterH$scores, N_CLUSTERS, individualFiles=T, "components")
	boxplotsMatrix (clusterH$properties, N_CLUSTERS, individualFiles=F, "properties")
	boxplotsMatrix (clusterH$properties, N_CLUSTERS, individualFiles=T, "properties")
	boxplotsMatrixAll (clusterH$scores)

	#boxplotsComponentsOnGroups (clusterH)
}

###################################################################################
## Preliminary PCA to filter data. First, full measures, second without the four distant
###################################################################################
preliminaryPCA <- function (values){
	pc = PC (values, graph=T, filename="pca-analysis-pre")
	sink ("pca-analysis-pre.output")
		print (pc$vectors, digits=2)
		print (pc$importance, digits=2)
	sink()

	#pca0 = principal(values, nfactors=N_COMPONENTS, residuals=F, rotate="none", scores=T, covar=T)
	pca0 = principal(values, nfactors=16, rotate="varimax")
	pca0$scores = as.matrix (values) %*% pca0$loadings
	outFilename = sprintf ("pca-analysis-pre-none-principal.output")
	sink (file=outFilename)
		print (pca0)
	sink()
	return (pc)
}
###################################################################################
# Full analysis process of the measures
###################################################################################
measureAnalysis <- function (valuesTraining, valuesTesting, clusterFilename, nComponents, label) {
	# preliminaryFirst16PCA  <<- getRotatedPrincipalComponents (allValues, 16, "-pca16cAllp")
	pca <- getRotatedPrincipalComponents (valuesTraining, nComponents, name("-pcaFinal",nComponents))
	clusterH <- makeClustering (clusterFilename, pca, valuesTraining)
	classification = classificateWholePaths (valuesTesting, clusterH$loadings, clusterH$centers, label)

	write.table (file="pca-cluster-centers.values", clusterH$centers, sep=" ", quote=F)
	write.table (file="pca-cluster-loadings.values", clusterH$loadings, sep=" ", quote=F)
	write.table (file="pca-classification.values", classification, sep=" ", quote=F)

	return (list (pca=pca, cluster=clusterH, classification=classification))
}

#-----------------------------------------------------------------------
# Do a Principal Components with rotated componets
# Return a PCA object with info of the pca analysis
# It perform PCA with three types of rotations: none, varimax and oblimin
#-----------------------------------------------------------------------
invisible (library(psych))  ## principal  -- PCA with rotated components
getRotatedPrincipalComponents <- function (values, N_COMPONENTS, label) {
	log ("Making PCA rotated NONE...")
	pca0 = principal(values, nfactors=N_COMPONENTS, residuals=F, rotate="none", scores=T, covar=T)
	pca0$scores = as.matrix (values) %*% pca0$loadings
	outFilename = sprintf ("pca-analysis%s-none.output", label)
	sink (file=outFilename)
		print (pca0)
	sink()

	log ("Making PCA rotated with varimax...")
	pca1 <<- principal(values,nfactors=N_COMPONENTS,residuals=F,rotate="varimax",scores=T,covar=T)
	pca1$scores = as.matrix (values) %*% pca1$loadings
	outFilename = sprintf ("pca-analysis%s-varimax.output", label)
	sink (file=outFilename)
		print (pca1)
	sink()

	## For check Sum of squares with none ratation
	log ("Making PCA rotated with oblimin")
	pca2 <- principal(values, nfactors=N_COMPONENTS, residuals=F, rotate = "oblimin", covar=T)
	outFilename = sprintf ("pca-analysis%s-oblimin.output", label)
	sink (file=outFilename)
		print (pca2)
	sink()

	return (list (loadings=pca1$loadings, scores=pca1$scores, pca0=pca0, pca1=pca1, pca2=pca2))
}

##############################################################################
# Print a log message with the parameter
###############################################################################
log <- function (...) {
	messages = unlist (list (...))
	cat (">>>>", messages, "\n")
}

###################################################################################
# Paste a string vector by collapsing the strings
###################################################################################
name <- function (...) {
	words = list(...);
	return (paste (unlist (words), collapse=""))
}

##############################################################################
# Create the output dir
# If it exists then rename it as old-XXXX
##############################################################################
createOutputDir <- function (outputDir) {
	moveToOldDir <- function (dir) {
		oldDir = name ("old-", dir)
		if (file.exists (oldDir)) 
			moveToOldDir (oldDir)
		cmm = sprintf ("mv %s %s", dir, oldDir)
		system (cmm)
	}

	if (file.exists (outputDir)) 
		moveToOldDir (outputDir)

	dir.create (outputDir)
}

##############################################################################
# Get names of the groups according to num of clusters
##############################################################################
getGroupNames <- function (N_CLUSTERS) {
	groupNames = paste (rep("g", N_CLUSTERS), 1:N_CLUSTERS, sep="")  # g1, g2, g3,   gK
	return (groupNames)
}
##############################################################################
# Get covariance
##############################################################################
getCovariances <- function (values, outputDir) {
	covariances = cov (values)
	write.table (round (covariances, digits=2), file="covariances.values", col.names=NA)
	return (covariances)
}

###########################################################
# Print and plot the MDS results
###########################################################
plotMDS <- function (mds, title, correlations=NULL, significances=NULL, dissimilarities=NULL) {
	outFilename = sprintf ("%s.pdf", title)
	pdf (outFilename, width=4, height=4)
		par (mar=c(4,4,2,2))
		x=mds[,1]; y=mds[,2]
		minX=min(x);maxX=max(x);minY=min(y);maxY=max(y)
		plot (x,y, type="n")
		text (x,y, rownames (mds))#, cex=1.2)
		lines (c(minX,maxX),c(0,0));lines (c(0,0),c(minY,maxY))
	dev.off()
}

###########################################################
# Run MDS with the simmilarities (non dissimilarites (default))
###########################################################
getMDS <- function (values, title) {
	#library(Hmisc)			# For significances (rcorr function)
	covariances = cov (values)
	#significances = rcorr (as.matrix (data), type="pearson")
	#dissimilarities = as.dist (1 - abs (correlations))
	#mds = cmdscale (dissimilarities, k=2)
	#plotMDS (mds, title, correlations, significances, dissimilarities)
	#similarities = dist (covariances)
	dissimilarities = as.dist (1 - abs (covariances))
	mds = cmdscale (dissimilarities, k=2)
	plotMDS (mds, title)
}
###########################################################
# MDS analysis
###########################################################
MDSAnalysis <- function (values, badProperties, label) {
	log ("MDS Analysis...")
	getMDS (values, name ("MDS-all", label))
	title = name ("MDS-no")
	for (prop in badProperties) {
		index = which (names (values) == prop)
		title = name (title, prop)
		values = subset (values, select = -index)
		badProperties = badProperties [badProperties != prop]
		getMDS (values, name (title, label))
	}
}

###############################################################################
# Return the order of groups according to folding levels.
# Currently it checks values of the NC and get he order from it
###############################################################################
getGroupOrder <- function (values, classes, N_CLUSTERS) {
	write.table (classes, file="classes-before-rename.values")
	nc = cbind (NC=values [,"NC"], classes)
	#nc = values [, c("NC","classes")]
	bx = boxplot (NC~classes, nc, plot=F)
	stats = bx$stats
	means = stats [3,]
	groupOrder = order (means)

	indices = list(N_CLUSTERS)
	n=1
	for (i in groupOrder) {
		indices [[n]] =  which (classes==i)
		n=n+1
	}

	for (i in 1:N_CLUSTERS)
		classes = replace (classes, indices [[i]], i)

	write.table (classes, file="classes-after-rename.values")
	return (classes)
}

###############################################################################
# Make clustering from scores of Principal Component Analysis
# Inputs: data values, number of Comps, number of Groups, link method, distance method
###############################################################################

clusteringHierarchical <- function (pca, values, savedData=NULL) {
	invisible (library (fpc))   # cluster.stats
	invisible (library (amap))  # hcluster
	invisible (library (cluster))  # hcluster
	distanceMeth = CLUSTER_DISTANCE
	linkMeth     = CLUSTER_LINKAGE
	scores = pca$scores 
	distances = dist (scores, method = distanceMeth) # Use also (best) "Dist" from (amap)

	# Make or load a Hierarchical Clustering and cut tree into 4 clusters
	if (is.null (savedData)) {
		# Another function which gives more info 
		# hcAgnes <<- agnes (scores, method="complete") 
		hc <- hcluster(scores, link=linkMeth, method=distanceMeth, nbproc=1) 
		classes <- cutree(hc, k=N_CLUSTERS)
		write.table (classes, file="classes-cuttree.values")
		stats = cluster.stats (distances, classes)
	}else{
		hc <-savedData$cluster
		classes <<- savedData$classes
		stats <- savedData$stats
	}
		
	# Rename classes according to heuristic criteria (NC order)
	#orderedClasses <-getGroupOrder (values, classes, N_CLUSTERS)
	orderedClasses <- classes
	# Get quality Silhoette info
	silhInfo = silhouette (orderedClasses, distances, do.clus.stat=FALSE, do.n.k=FALSE)
	# Add classes to scores and values
	scores = cbind (scores, classes=orderedClasses)	
	properties = cbind (values, classes=orderedClasses)
	#### Calculate score vector lengths (square root of the sum of squares of components)
	#### lengthsVct <- sqrt (apply (scores [,1:N_COMPONENTS], 1, function (x){sum(x^2)}))
	#### lengths <- cbind (length=lengthsVct, class=scores [,"classes"])
	# Calculate centers
	centers <- NULL
	for(k in 1:N_CLUSTERS)
		centers <- rbind(centers, colMeans(pca$scores[orderedClasses == k, , drop = FALSE]))

	log ("CENTERS:", round (centers, digits=3))
	
	return (list (scores=scores, loadings=pca$loadings, centers=centers, classes=orderedClasses, 
			cluster=hc, properties=properties, silhInfo=silhInfo, cluster=hc, stats=stats))
}

##############################################################################
# Make (or load) the hierarchical clustering using a complete linkage and euclidean distance
##############################################################################
makeClustering <- function (clusterFilename, pca, valuesTraining) {
	linkMeth = CLUSTER_LINKAGE
	distMeth = CLUSTER_DISTANCE
	log ("Clustering...", clusterFilename)
	# Make (or load) clustering using principal components
	clusterH = NULL
	if (file.exists (paste ("/tmp/", clusterFilename, sep=""))) {
		log ("Loading data of clustering ", clusterFilename)
		load (paste ("/tmp/", clusterFilename, sep=""))  # loading "clusterH" object
		clusterH <- clusteringHierarchical (pca, valuesTraining, clusterH)
	} else {
		log ("Creating data of clustering ", clusterFilename)
		clusterH <- clusteringHierarchical (pca, valuesTraining)
		log ("Saving data of clustering...")
		save ("clusterH", file=clusterFilename)
	}
	return (clusterH)
}

###############################################################################
# Plot in a pdf format 
###############################################################################
plotPdf <- function (plotCall, filename) {
	pdf (filename, title=filename, width=5, height=3)
	par(cex=0.8, xpd=T, mar=par()$mar+c(0,0,0,10))
	plotCall
	groupNames = c("g1: n=1020, s=.59", "g2: n=1175, s=.50", "g3: n=1662, s=.43", "g4: n=1542, s=0.60")
	colors = unique(1:N_CLUSTERS)
	legend ("right", legend=groupNames, col=unique (colors), inset=c(-0.7,0),
	        lty=1, lwd=1.3, title="Groups",cex=1.0, pch=1:N_CLUSTERS)
	dev.off()
}

###############################################################################
# Plot cluster and Silhoutte info for HIERARCHICAL cluters
###############################################################################
invisible (library(scatterplot3d)) #

plotClusteringInfo <- function (clusterInfo, N_CLUSTERS, label) {
	cluster = clusterInfo$cluster
	stats = clusterInfo$stats
	silhInfo = clusterInfo$silhInfo
	scores = clusterInfo$scores
	classes = clusterInfo$classes
	
	filename = name ("plotclust-", "3d", label, ".pdf")

	# Plot cluster in 3D
	filename = name ("plotclust-", "3d", label, ".pdf")
	pdf (filename, title=filename, width=4, height=4)
		par(cex=0.3, xpd=T)  # no legends
		scatterplot3d (scores[,1:3], color=classes, label.tick=F, mar=c(5,4,4,2)+0.2,
		               xlim=c(-0.8,0.8), zlim=c(-0.6,0.6), ylim=c(-0.8,0.8),
		               xlab="Stability", ylab="Compactness", zlab="Native-likeness")
		#legend ("topright", legend=groupNames, col=unique (colors), inset=c(-0.3,0),
	    #	    lty=1, lwd=1.5, title="Groups",cex=1.2, pch=1:N_CLUSTERS)
	dev.off()

	# Plot dendogram
	filename = name ("plotclust-", "tree", label, ".pdf")
	pdf (filename, title=filename, width=4, height=4)
		#par(cex=0.8, xpd=T, mar=par()$mar+c(0,0,0,10)) # with legends
		par(cex=0.7, xpd=T, mar=par()$mar+c(0,0,0,0))  # no legends
		plot(cluster, labels=F, main=NULL)
		classes = rect.hclust(cluster, k=N_CLUSTERS, border=c(1,4,2,3)) 

		filename = name ("plotclust-", "datapoings", label, ".output")
	dev.off()

	# Plot silhouette
	filename = name ("plotclust-", "silh", label, ".pdf")
	plotPdf (plot (silhInfo, col=1:N_CLUSTERS, main="Silhouete Info", do.clus.stat=F), filename)

	# Write statistics values of cluster (diameter, distance, separation, silhinfo, ...)
	filename = name ("plotclust", "stats", label, ".values")
	writeLines (paste (names(stats), stats, sep="="), filename)
}

###########################################################################
## Put the text of the five values next to the box
###########################################################################
textBoxplots <- function (bx,nBoxes, type="median") {
	stats = bx$stats
	statsNames <- c("mx","tp","md","lw","mn")

	points = c()
	for (i in 1:nBoxes) {
		
		sttMtx = matrix (stats [,i], 5,1)
		sttText = paste (statsNames, rev(round (sttMtx,3)), collapse="\n")
		sttMeans = sttMtx[3,]

		# Only media
	    statsNames <- c("")  #  No names
		sttMtx = matrix (stats [c(3),i], 1,1)  # Only the third value (media)
		points = append (points, stats [c(3),i])

		# Only min and max
	    #statsNames <- c("","")
		#sttMtx = matrix (stats [c(1,5),i], 2,1)
		#sttText = paste (statsNames, rev(round (sttMtx,3)), collapse="\n")
		# text (x=seq (1.5,nBoxes+1,by=1)[i], y=sttMeans, sttText, cex=0.7)

		# All labels and values
		#sttText = paste (statsNames, rev(round (sttMtx,3)), collapse="\n")
		#text (x=seq (1.5,nBoxes+1,by=1)[i], y=sttMeans, sttText, cex=0.7)

		# Only mean values
		#sttText = round (sttMeans, digits=2)

		sttText = paste (statsNames, rev(round (sttMtx,2)), collapse="\n")
		text (x=seq (1.4,nBoxes+1,by=1)[i], y=sttMeans, sttText, cex=0.8)
	}
	boxplotSmoothCurveLine (seq(1:3), points)
	#log ("Points:", points)
	#y = points
	#x = seq (1:3)
	#lo = loess (y~x)
	#xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
	#lines(xl, predict(lo,xl), col='red', lwd=2)

	##smoothingSpline = smooth.spline(x, y, spar=0.15)
	##plot(x,y)
	##lines(smoothingSpline)

	#plot (x,y)
	#lines(predict(lo), col='red', lwd=2)

}
textBoxplotsICF <- function (bx,nBoxes, type="median") {
	stats = bx$stats
	for (i in 1:nBoxes) {
		sttMtx = matrix (stats [,i], 5,1)
		sttMeans = sttMtx[3,]

		sttValues = round (sttMtx[c(3),1], 2)
		text (x=seq (1.4,nBoxes+1,by=1)[i], y=sttValues, sttValues, cex=1.0)
	}
}
###############################################################################
# Make boxplots and write their statistics for components split into groups
###############################################################################
boxplotsMatrix <- function (dataMatrix, N_CLUSTERS, individualFiles=T, typeOfData, label="") {
	filename = name ("boxplot-", typeOfData, label, ".pdf")
	groupNames = getGroupNames (N_CLUSTERS)  # g1, g2, g3,   gK

	if (!individualFiles) {
		nrows = ceiling ( (ncol (dataMatrix)-1)/3)
		pdf (filename, title=typeOfData, width=8.5, height=nrows*2.5); 
		par (mfrow=c(nrows, 3))
	}
	nBoxplots = ncol (dataMatrix) -1
	for (n in 1:nBoxplots) {
		boxTitle = colnames (dataMatrix) [n]
		filename = name ("boxplot-", typeOfData, "-", boxTitle, label, ".pdf")
		if (individualFiles) pdf (filename, title=boxTitle, width=4, height=5)

		bx = boxplot (dataMatrix[,n]~classes, dataMatrix, names=groupNames, boxwex=0.3, ylim=c(-0.7, 1.0))
		title (main=boxTitle, xlab = "Groups", ylab = "scores") 
		stats = bx$stats

		#textBoxplots (bx, N_CLUSTERS)

		if (individualFiles) dev.off()
	}
	if (!individualFiles) dev.off()
}
#------------------------------------------------------------------------------
# One plot with all boxplots of the groups on the components
#------------------------------------------------------------------------------
boxplotsMatrixAll <- function (dataMatrix) {
	groupNames = c ("g1","g2","g3","g4")

	filename = name ("boxplot-components-All.pdf")
	pdf (filename, title="plot", width=16, height=5)

	layout (t(1:3))
	par(oma=c(3.5, 3.5, 0.5, 0.5), mar=c(1,1,1,0), cex=1)

	nBoxplots = ncol (dataMatrix) -1
	for (c in 1:nBoxplots) {
		boxplot (dataMatrix[,c]~classes, dataMatrix, boxwex=0.3, ylim=c(-1,1), axes=F)

		axis (side=1, at=seq(4), labels=groupNames)
        Axis (side=2, labels=F)
		mtext (paste ("Component ",c), 3, 0, cex=1.5)
        if (c==1) {
       		axis(2, las=1)
       		mtext("Scores", 2, 3, cex=1.7)
		}
	}
	mtext ("Groups ", 1,2, outer=T, cex=1.7)
	dev.off()
}

#--------------------------------------------------------------------
# Plot a smooth line over xy points
#--------------------------------------------------------------------
boxplotSmoothCurveLine <- function (x,y) {
	lo = loess (y~x)
	xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
	lines(xl, predict(lo,xl), col='red', lwd=1)
}
#####################################################################
## Used by the boxplotsComponentsOnGroups 
## Label for multiple components and multiple groups
#####################################################################
gpcs = function (groupNames, namesPCs, scoresDF) {
	## Label observations of group "g" with its component "pc"
	lab = function (g, pc) {
		sb=subset (scoresDF, classes==g, select=pc);colnames (sb)=g;data.frame (sb, PC=pc)
	};
	## Label for multiple components 
	mpcs = function (g, namesPCs) {
		dfs=data.frame();for (p in namesPCs) {df=lab (g, p);dfs=rbind (dfs, df)};dfs 
	};

	lgp = list(); 
	for (g in groupNames){ 
		gpc = sapply(g, mpcs, namesPCs); 
		lgp [[g]] = data.frame (g=gpc[[1]], pc=gpc[[2]])
	};
	lgp
}
#####################################################################
## Create boxplot for the distribution of components on each group
## INPUT: Cluster data, and names for axis of PCs
#####################################################################
boxplotsComponentsOnGroups <- function (clusterH, namesPCs, namesGroups, 
                                        namesAxis, labelFile) {
	log (">>>> Creating boxplots for components on groups")
	#groupNames = getGroupNames (N_CLUSTERS)
	groupNames = c("g1", "g4", "g2", "g3")

	sc <- clusterH$scores
	sc [,2] = -1 * sc [,2]
	colnames (sc) = c (namesPCs, "classes")
	scoresDF = as.data.frame (sc)
	# Add "g" to the number of the level
	scoresDF[, "classes"]=paste("g",scoresDF [,"classes"], sep="") 

	lspcs <- gpcs (groupNames, namesPCs, scoresDF)

	# Create individual boxplots of distribution of components on groups
	k=1
	for (g in groupNames) {
		filename = name ("boxplot-groups-", labelFile, g,".pdf")
		log (filename)
		pdf (filename, title=g, width=4, height=4)
		par (cex=0.60)
		#tiff (filename, title=g, width=3, height=3, res=300, unit="in")
		#bitmap(filename, height = 2, width = 2, units = 'in', type="tifflzw", res=300)
		#tiff (filename, title=g, res=300)
			par(xpd=TRUE)
			#par(mar=c(4, 4, 2, 1 )) 
			#bx = boxplot (g~pc, lspcs[[g]], boxwex=0.3, ylim=c(-1,1))
			bx = boxplot (g~pc, lspcs[[g]], boxwex=0.3, ylim=c(-1,1))
			title (main=namesGroups[k], xlab = "Features", ylab = "Scores") 
			stats = bx$stats
			#log ("STATS:>>>>", dim(stats), stats)
			#textBoxplots (bx, length(namesPCs))
			#print (par()$mar)

			boxplotSmoothCurveLine (seq (1:3), bx$stat[3,])
			k = k + 1

		dev.off()
	}

	#
	# Make a single plot of all boxplots
	#
	filename = name ("boxplot-groups-", labelFile, "all.pdf")
	pdf (filename, title="plot", width=6.5, height=3.0)
		layout(t(1:4)) # 
		par(oma=c(3.5, 3.5, 0.5, 0.5), mar=c(6,1,1,0), cex=0.7)

		k=1
		for (g in groupNames) {
			bx = boxplot (g~pc, lspcs[[g]], boxwex=0.3, ylim=c(-1,1), xlim=c(0.7,3.7), axes=F)
			#axis (side=1, at=seq(3),labels=namesPCs, cex.axis=1.3)
			axis (1, at=seq(3),labels=F)
			text (x=seq(1:3), labels=namesPCs, par("usr")[3] - 0.4, srt=55, xpd=T, pos=1, cex=0.8)
			Axis(side=2, labels=F, cex=0.8)
			#mtext (paste ("Group ",k), 3,0, cex=1.5)
			#mtext (namesGroups [k], 3,0, cex=0.8)  # top 

			textBoxplots (bx, 3)
			if(g == "g1"){
				axis(2, las=1)
				mtext(namesAxis[2], 2, 3, cex=0.8)  # Left side
			}
			#if(g == "g1"){
			#	axis(2, las=1)
			#	mtext("Scores", 2, 3, cex=1.7)
			#	textBoxplots (bx, 3)
			#}
			k = k+1
		}
		mtext (namesAxis [1], 1,-1, outer=T, cex=0.8) # Botton line
	dev.off()

	# Make a boxplot of ICF Score on Folding Levels
	icfscoreCalc <- function (f1,f2,f3){ 
		icfs = (0.63*f1-0.09*f2+0.06*f3)/0.78;return (icfs) 
	}
	icfscore = mapply (icfscoreCalc, sc[,1],sc[,2],sc[,3]) 
	scICF = cbind (sc, ICF=icfscore)
	
	log ("Boxplots ICF")
	filename = name ("boxplot-groups-", labelFile, "ICF.pdf")
	pdf (filename, title="plot", width=8, height=5)
		bx = boxplot(ICF~classes, scICF, names=namesGroups, col=c(8,2,3,4),
		     xlab="Folding Levels", ylab="Score", main="ICF Score", boxwex=0.4)
		textBoxplotsICF (bx, 4, "ICF")
	dev.off()
	filename = name ("boxplot-groups-", labelFile, "ICF-stats.values")
	write.table (file=filename, bx$stats)
}
###############################################################################
# Analysis of Components
# Make boxplots of the different PCAs made in cross validation
###############################################################################
boxplotsCrossValPCAs <- function (label) {
	statsNames <- c("max","top","med","low","min")
	pcas = read.table ("inputs/pcas-loadings-crossval-20p-5g.values", header=T)

	allPCs = list (pc1=c("NC","HB","CS","AS","DF"), pc2=c("RG","RM","LR","SS"), pc3=c("CO","SA"))
	n = 0
	for (pc in allPCs) {
		n = n + 1
		filename = name ("boxplot-crossval", label, ".pdf")
		pdf (filename, title=filename,  width=5.5, height=4.25)
		nK = length (pc)

		props3 = pcas [pcas$Pr %in% pc,]
		props3 = cbind(abs (props3[,1:3]), Pr=as.vector(props3[,4]))
		bx = boxplot (props3 [, n]~Pr, props3, main = sprintf ("PC%s",n), boxwex=0.3)
		for (i in 1:nK) {
				sttMtx = matrix (bx$stats [,i], 5,1)
				sttText = paste (statsNames, rev(round (sttMtx,3)), collapse="\n")
				sttMeans = sttMtx[2,]  
				text (x=seq (1.4,nK+1,by=1)[i], y=sttMeans, sttText, cex=0.6)
		}
		dev.off()
	}
}
###############################################################################
# Classificate a specific pathway using clustering centers and PCA loadings
# Cluster calculation for pathway conformations
# Classification of one pathway into the valuesPath 
###############################################################################
classificateOnePath <- function (valuesPath, loadings, centers) {
	#log ("Classificating...")
	# Calculating scores, transforming data to zero mean and unit variance
	valuesPath.std <- as.matrix (valuesPath)

	# Calculate scores from eigen analysis
	valuesPath.scores <- valuesPath.std %*% loadings
	valuesPath.scores <- valuesPath.scores [,1:N_COMPONENTS]

	# Assign classes comparing centers with each score
	classes <- vector (length= nrow (valuesPath.scores))
	for (i in 1:nrow (valuesPath.scores)) {
		sc <- valuesPath.scores [i,]
		centersSc <- rbind (centers, sc)
		matrixDistances <- as.matrix (dist (centersSc, method=CLUSTER_DISTANCE))
		distances <- matrixDistances ["sc",][1:N_CLUSTERS] ## 4 groups of the clustering
		classes [i] <- which.min (distances)
	}
	# Calculate score vector lengths (square root of the sum of squares of components)
	lengths <- sqrt (apply (valuesPath.scores, 1, function (x){sum(x^2)}))
	
	# Create a big matrix with previous results
	classification <- cbind (valuesPath.scores, classes, lengths)
	componentNames =  paste ("c",seq (1:N_COMPONENTS), sep="")
	colnames (classification) <- c (componentNames,"classes","lengths")
	return (data.frame (classification))

}
###############################################################################
# Plot pathway: components and groups
# Returns the matrix of classification of the test "values"
# Classification whole pathways into the values 
###############################################################################
classificateWholePaths <- function (values, loadings, centers, label) {
	log ("Classificating and Plotting Pathways...")
	pathwayNames = unique (unlist (Map (function (i) i[1],strsplit (rownames (values), "-"))))
	log (">>> ", "Evaluating", pathwayNames)
	allPathsClassification = data.frame ()
	for (path in pathwayNames) {
		index = grep (path, rownames (values))
		valuesPath = values  [index,]

		# Classificate the pathway using clustering centers and PCA loadings
		classification <<- classificateOnePath (valuesPath, loadings, centers)

		scoresPath = classification [,1:N_COMPONENTS]
		#classesPath = classification [,"classes", drop=F]
		classesPath = classification [,"classes"]
		lengthsPath = classification [, "lengths"]

		scoresWithGroups = cbind (scoresPath, group=classesPath)
		write.table (scoresWithGroups, file=name ("scores-test-", path, ".values"),quote=F)

		# Added to get a table with all values and scores
		outFilename = name ("pathway-", path,".values")
		write.table (round (valuesPath,3), file=outFilename, quote=F)

		plotPathway (N_CLUSTERS, classesPath, path, "")
		
		allPathsClassification = rbind (allPathsClassification, classesPath)
	}
	write.table (file="all-paths-classification.values", allPathsClassification)
	return (allPathsClassification)
}

###############################################################################
## Plot pathway with its conformations showing its groups (optional)
###############################################################################
plotPathway <- function (N_CLUSTERS, classesPath, pathwayName, label) {
	filename =  name ("plotpath-", pathwayName, label, ".pdf")
	pdf (filename, title=filename, width=5, height=3)
		#par(mar = c(4.0,0.01,2.5,0.01))
		par(mar = c(4.0,4.01,2.5,0.01), mgp=c(3.0,1,0))
		LWD=1.5; CEX=0.7
		N = length (classesPath)
		time = 1:N	# Time to print on x axis
		title = pathwayName
		#groupNames = c ("g1","g2","g3","g4")
		#groupNames = c ("U:Unfolded","E:Early Intermediate","L:Late Intermediate","F:Folded")
		groupNames = c ("Unfolded","Early\nInterm.", "Late\nInterm.","Folded")
		strClassesPath = paste ("g", classesPath, sep="")

		plot (x=time, y=classesPath, col="black", xlab="Path Step", 
			main=paste ("Predictions for", pathwayName), 
			ylab="Folding Levels", lty=1, type="l", cex=CEX, axes=F, cex.axis=1.0)

		par (cex = CEX)
		axis (1, at=seq (0,N+10, 10))
		axis (2, at=1:N_CLUSTERS, lab=groupNames, las=2)
		#axis (2, at=1:N_CLUSTERS, lab=1:N_CLUSTERS)
		#axis (3, at=seq (1,N,1), lab=strClassesPath, cex.axis=0.8, line=-1, tick=F)

		for (i in seq (1,N+1,2)) 
				points (x=time[i], y=classesPath[i], col=classesPath[i], 
		            lty=1, lwd=LWD, type="o", pch=c(classesPath[i]),cex=CEX)

		#text (x=15, y=3, label="Path Step", pos=1, cex=0.5)

		#legend ("bottomright", legend=groupNames, col=unique (1:N_CLUSTERS), 
		#        lty=1, lwd=LWD, title="Folding Levels",cex=1.2, pch=1:N_CLUSTERS)
	dev.off()
}

###################################################################################
# Removes properties not included in the analysis
###################################################################################
cleanValues  <- function (values, BAD_PROPERTIES) {
	PROPS = !(names (values) %in% BAD_PROPERTIES) 
	values = values [PROPS]
	return (values)
}
#############################################################################
## Split the dataset "x" of names in "n" subsets (training and test) 
## Return a list with "n" training and testing subsets
#############################################################################
getTrainingTestingDatasets <- function(x,n) {
	tst = split(x, factor(sort(rank(x)%%n)))
	trn = sapply (tst, function (y) { x [ ! x %in% unlist(y)] }, simplify=F)
	return (list (training=trn, test=tst))
}

###################################################################################
# Return the full values corresponding to "names" 
###################################################################################
getValuesDataset  <- function (values, names) {
	log ("Getting values dataset...")
	# get the full values given partial row name
	filter <- function (name) {values [name,]}

	datasetLists = sapply (names, filter, USE.NAMES=T, simplify=F)
	df = data.frame ()
	datasetDataframe = sapply (datasetLists, function (lst) {df <<- rbind (df, lst)})
	
	return (df)
}
getValuesPathways  <- function (values, names) {
	# get the full values given partial row name
	filter <- function (path) {values [grep (path[[1]], rownames (values)),]}

	datasetLists = sapply (names, filter, USE.NAMES=T, simplify=F)
	df = data.frame ()
	datasetDataframe = sapply (datasetLists, function (lst) {df <<- rbind (df, lst)})
	
	return (df)
}

###################################################################################
# Split the dataset in two sets :training and testing, for cross validation
# The number of testing observations is given using the "dilutionFactor" (0 < df < 1  )
# Return a list with the training and testing sublist for each subdataset
#################################################measures##################################
splitDataset  <- function (dataset, dilutionFactor, numberOfSets=NULL) {
	## Define chunks parameters
	sizeDataset    = length (dataset)
	sizeTraining   = floor (dilutionFactor * sizeDataset)
	sizeTesting    = sizeDataset - sizeTraining
	numberOfChunks = floor (sizeDataset / sizeTesting)

	testing = split(dataset, factor(sort(rank(dataset)%%numberOfChunks)))
	training = sapply (testing, function (y) {dataset[! dataset %in% unlist(y)] }, simplify=F)

	log ("SPLIT Size",sizeDataset,"Factor",dilutionFactor,"N-chunks", numberOfChunks)
	log ("Size-test:" , sizeDataset/numberOfChunks,"Size-train:" , sizeDataset - sizeDataset/numberOfChunks)

	if (!is.null (numberOfSets)) {
		testing = testing [[1]]
		training = training [[1]]
	}
	#print (str (testing))
	#print (str (training))

	return (list (trainingNames=training, testingNames=testing))
}

###################################################################################
# Calculate the percentage of coincidence over the k-fold cross validation
# 100% full, 75% more coinc than diffs, 50% 1/2, and 25% 
# The result is added as a column to the input matrix
###################################################################################
getStimationsForPredictions <- function (classes) {
	# Round % coincidence to thesholds (25%, 50%, 75%, 100%)
	qualityLabel <- function (stimation) {
		if (stimation < 0.25)      return (25)
		else if (stimation < 0.50) return (50)
		else if (stimation < 0.75) return (75)
		else                       return (100)
	}
	nFolds = ncol (classes)
	nClasses = N_CLUSTERS
	stimationMatrix = data.frame ()
	
	for (i in 1:nrow (classes)) { 
		# Count of agreements, according to folding level, for the k-fold rounds
		stimationVector = mapply (function (j) length (which (classes [i,]==j)), 1:nClasses)
		# Stimate the rate for each folding level agreement
		stimationRate =  mapply (function (k) stimationVector [k] / nFolds, 1:nClasses)
		# Select the maximal stimation agreement
		stimationVector ["Stimation"] = round (max (stimationRate), digits=2)
		stimationVector ["Quality"] = qualityLabel (max (stimationRate))

		# Builds a matrix of stimations
		stimationMatrix = rbind (stimationMatrix,  stimationVector)
	}
	# Label de stimation on each pathway's conformation
	colnames (stimationMatrix) = c("L1", "L2", "L3", "L4", "Stimation", "Quality")
	classes = cbind (classes,  Quality=stimationMatrix[,"Quality"])
	stimations=mapply (c(25,50,75,100), FUN=function (i) sum (classes [,"Quality"]==i))
	names (stimations) = c("D","C","B","A")

	return (stimations)
}

plotStimationsForPredictions <- function (stimations, classes) {
	# Pie Chart with Percentages
	slices <- subset (stimations, stimations>0)
	lbls <- names (slices)
	pct <- round(slices/sum(slices)*100)
	lbls <- paste(lbls, pct) # add percents to labels
	lbls <- paste(lbls,"%",sep="") # ad % to labels
	pdf ("stimations.pdf" )
	pie(slices,labels = lbls, col=rainbow(length(lbls)), main="Stimations") 
	dev.off()
}
###################################################################################
## Return a list with the observations corresponding to each group of clustering
###################################################################################
getNeighborsGroups <- function (classes) {
	listGroups =  sapply (1:N_CLUSTERS, function (x) subset (classes, classes==x), simplify=F)
	names (listGroups) = paste ("g", 1:N_CLUSTERS, sep="")
	return (listGroups)
}

###################################################################################
## Mark the neighbors for each vector corresponding to each group
###################################################################################
setNeighborsVector <- function (classification, testingNames, nCols) {
	log ("Setting neighbors in a vector...")
	classes = subset (classification, select="classes")
	groups = getNeighborsGroups (classes)
	
	## Cols for observations, rows for groups
	nRows  = length (groups) 

	mat = matrix (rep (0, nRows*nCols), nrow=nRows, ncol=nCols, dimnames=list(1:nRows,1:nCols))
	mat [,testingNames] = 1

	for (i in 1:nRows) {
		g = groups [[i]]
		posNames = as.integer (rownames (g))
		mat [i,posNames] = 1
	}

	return (mat)
}

###################################################################################
## Set the group level to the neighborgs of each group
###################################################################################
setNeighborsMatrix <- function (classification, testingNames, nRows=NULL) {
	log ("Setting neighbors in a matrix ...")
	classes = subset (classification, select="classes")
	groups <<- getNeighborsGroups (classes)

	if (is.null (nRows) )
		nRows = nrow (classes)

	m = length (groups)

	mat = matrix (rep (0, nRows*nRows), nrow=nRows, ncol=nRows, dimnames=list (1:nRows,1:nRows)) 
	mat [testingNames,] = 1
	mat [,testingNames] = 1
	for (i in 1:m) {
		g = groups [[i]]
		names = rownames (g)

		mat [names,names] = 1
	}
	return(mat)
}

###################################################################################
## Return the names of the pathways in a dataset of values of proteins
## For the dataset of Parasol, proteins name has a prefix the pathway name (e.g 2YCC-p01.pdb)
###################################################################################
getPathwayNames <- function (values) {
	pathwayNames = unique (unlist (Map (function (i) i[1],strsplit (rownames (values), "-"))))
	return (pathwayNames)
}

###################################################################################
## Evaluate the resampling rounds
## The inputs are the matrix or vector of clustering results for the 
## full dataset and the  ensemble of resamplings 
###################################################################################
evaluateResampling<- function (fullDataset, ensembleDataset, nSamples, dilutionFactor, label) {
	log ("Evaluation of resampling...")
	mfSum = sum (fullDataset)
	mRows = length (ensembleDataset)
	mArray = c()
	for (i in 1:mRows) {
		miMul = fullDataset * ensembleDataset [[i]]
		mRatio = sum (miMul) / mfSum
		mArray = append (mArray, mRatio)
		log (i, mRatio)
	}
	outputFilename = sprintf ("validation-n%s-df%s-%s.results", nSamples, dilutionFactor*100, label)
	sink (file=outputFilename)
	cat ("Method:", label, "\nN samples:", nSamples, "\nDilution Factor:", dilutionFactor)
	cat ("\nN rounds:", length(mArray), "\nMean resamples:", mArray, "\nTotal mean:", mean (mArray))
	sink()

	return (mArray)
}

###################################################################################
# Cross validation by resampling
# Cluster all the samples and compare with cluster of n resamplings
###################################################################################
measureAnalysisResamplingVersion <- function (allValues, nSamples, dilutionFactor) {
	filteredValues   = cleanValues (allValues, BAD_PROPERTIES)
	values           = filteredValues [sample(rownames(filteredValues), size=nSamples),]

	## Change names to protein samples to integer 1:N numbers
	rownames (values)   <- 1:nrow (values)

	## Full clustering (100%)
	clusterFilename = sprintf ("cluster%sSample.RData", nSamples)
	pca <- getRotatedPrincipalComponents (values, "")
	clusterH = newMakeClustering (clusterFilename, pca, values)
	classification = classificateOnePath (values, clusterH$loadings, clusterH$centers)

	nRows = nrow (classification)

	fullNeighborsMatrix  = setNeighborsMatrix (classification, c())
	fullNeighborsVector  = setNeighborsVector (classification, c(), nrow (classification))

	ensembleNeighborsMatrix = list()
	ensembleNeighborsVector = list()

	# Partial clustering (resampling of m subsamples)
	log ("\n>>>> Resampling\n")
	datasets = splitDataset (rownames (values), dilutionFactor) 
	for (i in 1:length (datasets$trainingNames)) {
		log (">>>> Resampling: ", i) 
		valuesTraining = getValuesDataset (values, datasets$trainingNames[[i]])
		testingNames  = datasets$testingNames[[i]]

		clusterFilename = sprintf ("cluster%s-%s-Sample.RData", nSamples,i)
		pca <- getRotatedPrincipalComponents (valuesTraining, "")
		clusterH = newMakeClustering (clusterFilename, pca, valuesTraining)
		classification = classificateOnePath (valuesTraining, clusterH$loadings, clusterH$centers)

		partialNeighborsMatrix   = setNeighborsMatrix (classification, testingNames, nrow (values))
		partialNeighborsVector   = setNeighborsVector (classification, testingNames, nrow (values))

		ensembleNeighborsMatrix [[as.character (i)]] = partialNeighborsMatrix
		ensembleNeighborsVector [[as.character (i)]] = partialNeighborsVector
	}

	evalMatrix = evaluateResampling(fullNeighborsMatrix, ensembleNeighborsMatrix, nSamples, dilutionFactor, "Matrix")
	evalVector = evaluateResampling(fullNeighborsVector, ensembleNeighborsVector, nSamples, dilutionFactor, "Vector")
}

#####################################################################
# Calculate the PCA given an input data
#####################################################################
PC<-function(X,method="eigen",scaled=T,graph=F,rm.na=T,print.results=T, filename=""){
	if (any(is.na(X))){
		tmp<-X
		if(rm.na==T){X<-na.omit(data.frame(X));X<-as.matrix(X)}
		else{X[is.na(X)] = matrix(apply(X, 2, mean, na.rm = TRUE),
			 ncol = ncol(X), nrow = nrow(X), byrow = TRUE)[is.na(X)]}}
	else{tmp<-X}

	if(method=="eigen"){
	if(scaled==1){X1=cor(X);X2=scale(X)}
	else{X1=cov(X);X2=scale(X,scale=F)}

	total.var<-sum(diag(cov(X2)))
	values<-eigen(X1)$values;vectors<-eigen(X1)$vectors;sdev=sqrt(values)}

	if(method=="svd"){
		if(sum(scaled,center)>1){X2<-scale(X)}
		else{if(scaled==1){X2=scale(X,center=F)}	
			 else{if(center==1){X2=scale(X,scale=F)}
				  else{X2=X}}}

	total.var<-sum(diag(cov(X2)))
	var<-nrow(X2)-1
	vectors<-svd(X2)$v;sdev=svd(X2)$d/sqrt(var);values<-sdev*sdev}
	prop.var<-rep(NA,ncol(X));cum.var<-rep(NA,ncol(X));scores<-X2%*%vectors
	namex<-as.character(1:ncol(X));scorenames<-rep(NA,ncol(X))

	for(i in 1:(ncol(X))){
		scorenames[i]<-do.call(paste,c("PC",as.list(namex[i]),sep=""))
		}
	colnames(scores)<-scorenames
	rownames(vectors)<-colnames(X);colnames(vectors)<-scorenames
	for(i in 1:ncol(X)){prop.var[i]<-var(scores[,i])/total.var}
	for(i in 1:ncol(X)){cum.var[i]<-sum(prop.var[1:i])}

	importance<-t(matrix(c(sdev,prop.var,cum.var),ncol=3))
	importance<-as.table(importance)
	colnames(importance)<-scorenames
	rownames(importance)<-c("Standard Deviation","Proportion of Variance","Cumulative
	Proportion")
	z<-list(values=values,vectors=vectors,scores=scores,importance=importance
	,sdev=sdev)
	if(graph==1){
		pdf (file=sprintf ("%s-biplot.pdf", filename))
			#par (mar=c(2.8,2.8,1,1), cex=0.7, mgp = c(1.5,0.8,0)) # Outer margins, scale fonts, pos labels (
			par (mar=c(4.5,4.5 ,2.5,2), cex.axis=1.7,cex.lab=1.7)
			xlabel = do.call (paste,c("C1 (",as.list(round(z$importance[2,1]*100,2)),"%)",sep=""))
			ylabel = do.call(paste,c("C2 (",as.list(round(z$importance[2,2]*100,2)),"%)",sep=""))
			biplot(scores[,1:2],vectors[,1:2], xlab=xlabel, ylab=ylabel, cex=c(0.15, 1.5), ylim=c(-5,4),xlim=c(-8,4))
		dev.off()

		pdf (file=sprintf ("%s-screeplot.pdf", filename))
			#par (mar=c(2.8,2.8,1,1), cex=0.7, mgp = c(1.5,0.8,0)) # Outer margins, scale fonts, pos labels (
			par (mar=c(4.5,4.5 ,2.5,2), cex.axis=1.7,cex.lab=1.7)
			screeplot(z,type='l', main=NULL)
			abline(1,0,col='red',lty=2)
			title (xlab="Components")
		dev.off()

		pdf (file=sprintf ("%s-groups.pdf", filename))
			sc = scores
			par (mar=c(4.5,4.5 ,2.5,2), cex.axis=1.7,cex.lab=1.7)
			plot (sc[,1],sc[,2], type="n", ylab="C2",xlab="C1")
			frac = nrow (sc)/3
			r1=1:frac; r2=(frac+1):(2*frac); r3=(2*frac+1):(3*frac)
			points (sc[r1,1],sc[r1,2], col="red", pch=1)
			points (sc[r2,1],sc[r2,2], col="blue", pch=2)
			points (sc[r3,1],sc[r3,2], col="green", pch=3)
		dev.off()
	}
	if(print.results==T){
		if(method=="eigen"){print("PCA Analysis Using Spectral Decomposition")}
		if(method=="svd"){print("PCA Analyis Using Singular Value Decomposition")}
		if (any(is.na(tmp))){
		if(rm.na==T){print("Warning:One or more rows of data were omitted from
		analysis")}
		if(rm.na==F){print("Warning: Mean of the variable was used for Missing
		values")}}
		print(importance)
	}
	z<-list(values=values,vectors=vectors,scores=scores,importance=importance
	,sdev=sdev)
}


###################################################################################
# MAIN 
# Run measure analysis by selecting various datasets for training and test (cross validattion)
###################################################################################
main ()


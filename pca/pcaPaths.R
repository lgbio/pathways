
options(width=500)

path="evaluaciones-normalizadas/"
filenames = list.files(path)
filenamesFullpath=paste (path, filenames, sep="")

# create a full dataset of all values of pathways 
nPathways = length (filenames)
filesDataframe = data.frame()
for (i in 1:nConformations) {
	df <- read.table (filenamesFullpath[i], row.names=1);
	filesDataframe <- rbind (filesDataframe, df [1:5,])
}
df = filesDataframe

crr=cor (filesDataframe)

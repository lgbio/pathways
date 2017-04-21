#!/usr/bin/python

USAGE  = "Calculate the structural similarity between structures using structural alphabets"
USAGE += "USAGE: structural.py <Target PDB >  <Reference PDB> [outputFilename]\n"

import os, sys
import uuid

TMPLABEL = "tmp_%s_" % str (uuid.uuid4())

#-------------------------------------------------------------
# MAIN 
#-------------------------------------------------------------
def main (args):
	# CHECK ARGUMENTS 
	if len (sys.argv) < 3:
		print USAGE
		sys.exit (0)
	elif len (sys.argv) == 4:
		sys.stdout = open (sys.argv[3], "w")

	# Load the alphabet
	alphabet = os.getenv ("EVAL_ALPHABET")

	# Process the target 
	pdbTarget = sys.argv [1]
	outputTarget = "tmp_%s_%s.output" % (TMPLABEL, name (pdbTarget))
	logTarget = "tmp_%s_%s.log" % (TMPLABEL, name (pdbTarget))

	params = {"IN": pdbTarget, "ALPHABET": alphabet, "OUT": outputTarget, "LOG": logTarget }
	cmm1 =  "encoder -g -F %(IN)s -A %(ALPHABET)s -T l -O %(OUT)s &> %(LOG)s" % params
	#print (cmm1)
	os.system (cmm1)

	# Process the reference 
	pdbRef = sys.argv [2]
	outputRef = "tmp_%s_%s.output" % (TMPLABEL, name (pdbRef))
	logRef = "tmp_%s_%s.log" % (TMPLABEL, name (pdbRef))

	params = {"IN": pdbRef, "ALPHABET": alphabet, "OUT": outputRef, "LOG": logRef }
	cmm2 =  "encoder -g -F %(IN)s -A %(ALPHABET)s -T l -O %(OUT)s &> %(LOG)s" % params
	#print (cmm2)
	os.system (cmm2)

	# Get the outputs
	encodeTarget = open (outputTarget).readlines()[2]
	encodeRef = open (outputRef).readlines()[2]

	# Evaluate the similarity
	equals = len (filter (lambda x:x[0]==x[1], zip (encodeTarget, encodeRef)))
	proportion = equals / float (len (encodeRef))

	# Clean temporary files
	clean ([outputTarget, outputRef, logTarget, logRef])
	print proportion

###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
	baseName = os.path.basename (namefile)
	newNamefile = baseName.split(".") [0]

	return newNamefile

##################################################################
# Delete temporal files created during the process
##################################################################
def clean (listOfFiles):
	log (listOfFiles)
	for file in listOfFiles:
		if os.path.lexists (file):
			os.system ("rm " + file)

#-------------------------------------------------------------
# Print a message to the error output stream
#-------------------------------------------------------------
def log (message):
	string=">>> StructSim: "
	if type (message) != list:
		string +=  str (message)
	else:
		for i in message: string += str (i) + " "
	
	sys.stderr.write (string+"\n")
#-------------------------------------------------------------
# Define the  output for errors and log messages
#-------------------------------------------------------------
def defineMessagesOutput ():
	stderr = os.getenv ("EVAL_STDERR")
	if stderr == None:
		sys.stderr = sys.stdout
	else:
		sys.stderr = open (stderr, "a")

#-------------------------------------------------------------
# Call to main
#-------------------------------------------------------------
if __name__ == "__main__":
	defineMessagesOutput ()
	main (sys.argv)


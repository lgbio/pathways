#!/usr/bin/python

import os, sys

args = sys.argv

oldPath = args [1]
newPath = args [2]

fileLst = os.listdir (".")

for file in fileLst:
	if os.path.islink (file):
		print file, 
		sl = os.readlink (file)

		newLink = sl.replace (oldPath, newPath)
		os.unlink (file)
		print ">", newLink
		os.symlink (newLink, file)


#!/usr/bin/python

import os, sys

args = sys.argv

if len (args) < 2:
	USAGE="%s <oldStr> <newStr>" % args[0]
	print USAGE
	sys.exit (0)

oldStr = args [1]
newStr = args [2]

filesLst = os.listdir (".")

for f in filesLst:
	if os.path.islink (f):
		oldLink = os.readlink (f)
		newLink = oldLink.replace (oldStr, newStr)
		
		print ">>>", f, oldLink, newLink
		os.unlink (f)
		os.symlink (newLink, f)

#!/usr/bin/python

# Split like python split but the input is redirected "<"

import sys

delimiter = sys.argv[1]
field = int (sys.argv[2])
string = raw_input()

"""
string = params.split(' "')[0]
print string
tmp = params.split('"')[1] 
delimiter = tmp.split('"')[0] 
field = params.split()[-1]
print field
"""

if delimiter == "":
	print string.split ()[int (field)]
else:
	print string.split (delimiter)[int (field)]

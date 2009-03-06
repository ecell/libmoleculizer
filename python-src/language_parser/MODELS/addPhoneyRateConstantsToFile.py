#!/usr/local/bin/python

import getopt, sys

PHONEYRATECONSTANT = 1e5

filename = ""
options, arguments = getopt.getopt(sys.argv[1:], "f:", ["file="])
for opt, atr in options:
    if opt in ("--file", "-f"):
	filename = atr
if filename == "":
    sys.exit(0)

lines = open(filename).readlines()
lines = [x.strip() for x in lines]

outputlines = []

for line in lines:
    if line != "" and line[-1] == "=":
	line = line + str(PHONEYRATECONSTANT)
    outputlines.append(line)


outputtext = '\n'.join(outputlines)
output = open(filename, 'w')
output.write(outputtext)
output.close()

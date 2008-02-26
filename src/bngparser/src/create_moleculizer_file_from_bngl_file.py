#!/usr/bin/python

###############################################################################
# BNGMZRConverter - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007 Nathan Addy
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Contact information:
#   Nathan Addy, Scientific Programmer	Voice: 510-981-8748
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#   2168 Shattuck Ave.                  
#   Berkeley, CA 94704
###############################################################################

import getopt, sys, re, pdb
import moleculizer

def main():
    # Handle command line options
    options = {}
    options["inputfile"] = ""
    options["outputfile"] = "default.xml"
    options["events_and_streams_mixin"] = False

    processOptions(options)
    if not ensureAppropriateOptions(options):
        usage()

    bngFile = open(options["inputfile"]).readlines()

    try:
        parameterBlock, moleculeTypeBlock, seedSpeciesBlock, reactionRuleBlock = \
                        parseBlockTypesFromBNG(bngFile)
    except Exception, e:
        barf("Error, bng file could not be parsed.")

    outputFile = moleculizer.MoleculizerObject(options["outputfile"])

    outputFile.addParameterBlock(parameterBlock)
    outputFile.addMoleculeTypeBlock(moleculeTypeBlock)
    outputFile.addSeedSpeciesBlock(seedSpeciesBlock)
    outputFile.addReactionRuleBlock(reactionRuleBlock)

    if options["events_and_streams_mixin"]:
        outputFile.addStreamsAndEventsMixin(options["events_and_streams_mixin"])

    outputFile.write()
    outputFile.close()

def processOptions(optionsDict):
    try:
	options, arguments = getopt.getopt(sys.argv[1:], "hf:o:v:m:", ["help", "file=", "output=", "version", "mixin="])
    except getopt.error, msg:
	print msg
	print "For help, use --help"
	sys.exit(0)

    for opt, atr in options:
	if opt in ("--help", "-h"):
	    usage()
            sys.exit(0)
	if opt in ("--file", "-f"):
	    optionsDict["inputfile"] = atr
	    fileRequired = True
	if opt in ("--output", "-o"):
	    optionsDict["outputfile"] = atr
	if opt in ("--version", "-v"):
	    version()
        if opt in ("--mixin", "-m"):
            optionsDict["events_and_streams_mixin"] = atr

def ensureAppropriateOptions(options):
    if not options["inputfile"]:
        return False
    else:
        return True

def usage():
    print "create_moleculizer_file_from_bngl_file.py version 0.5"
    print "\t-- Creates a Moleculizer input file from a BioNetGen input file"
    print
    print "Usage:"
    print "create_moleculizer_file_from_bngl_file.py -f inputfile [-o outputfile] [-m eventsandstreamsmixin]"
    print "create_moleculizer_file_from_bngl_file.py [-h] or [--help]"
    print 
    print "Nathan Addy (addy@molsci.org)"
    print "Molecular Sciences Institute, Berkeley, CA"
    print "June 22, 2007"
    print
    print 
    print "This is a preliminary version made for converting the bngl models produced through the alpha wiki."
    print "A major limitation of this version is that it requires a 'molecule types' block as well as a 'seed species'"
    print "block.  It then parses both 'species' and 'seed species' blocks in the same way.  I.e. implicit"
    print "molecule declaration is not allowed."

def version():
    print "BNGMZRConverter 0.5"
    print "Nathan Addy (addy@molsci.org)"
    print "Molecular Sciences Institute, Berkeley, CA"
    print "June 22, 2007"
    print 
    print "Licensed under GPL v2."
    print 

    print "This program is free software; you can redistribute it and/or modify"
    print "it under the terms of the GNU General Public License as published by"
    print "the Free Software Foundation; either version 2 of the License, or"
    print "(at your option) any later version."
    print
    print "This program is distributed in the hope that it will be useful,"
    print "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    print "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    print "GNU General Public License for more details."
    print
    print "You should have received a copy of the GNU General Public License"
    print "along with this program; if not, write to the Free Software"
    print "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA"
    print 
    print "Contact information:"
    print "Nathan Addy, Scientific Programmer"
    print "The Molecular Sciences Institute"
    print "2168 Shattuck Ave.                  Email: addy@molsci.org"
    print "Berkeley, CA 94704"

def parseBlockTypesFromBNG(bngFile):
    bngFile = [re.sub("#.*$", "", x) for x in bngFile] # Delete all comments
    bngFile = [x.strip() for x in bngFile]
    bngFile = [x for x in bngFile if x != ""] # This must be last, because line.strip() results in some empty lines.

    bngFile = '\n'.join(bngFile)
    bngFile = re.sub(r"\\\s*\n\s*", " ", bngFile)
    bngFile = bngFile.split("\n")

    try: 
        paramStartNdx = bngFile.index("begin parameters") + 1
        paramEndNdx = bngFile.index("end parameters")

    except ValueError,e:
        printerror("Could not locate parameter block")
        raise e

    try:
        molStartNdx = bngFile.index("begin molecule types") + 1
        molEndNdx = bngFile.index("end molecule types")
    except ValueError,e:
        printerror("Could not locate molecule block")
        raise e

    try:
        seedStartNdx = bngFile.index("begin seed species") + 1
        seedEndNdx  = bngFile.index("end seed species")
    except ValueError, e:
        printerror("Could not locate seed species block")
        raise e

    try:
        rxnStartNdx = bngFile.index("begin reaction rules") + 1
        rxnEndNdx = bngFile.index("end reaction rules")
    except ValueError, e:
        printerror("Could not locate reaction rule block")
        raise e

    parameterBlock = bngFile[paramStartNdx:paramEndNdx]
    moleculeBlock = bngFile[molStartNdx:molEndNdx]
    seedBlock = bngFile[seedStartNdx : seedEndNdx]
    rxnBlock = bngFile[rxnStartNdx:rxnEndNdx]

    return parameterBlock, moleculeBlock, seedBlock, rxnBlock

def barf(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.write("Crashing....\n")
    sys.exit(1)

def printerror(msg):
    sys.stderr.write(msg + '\n')
    return

if __name__ == "__main__":
    main()

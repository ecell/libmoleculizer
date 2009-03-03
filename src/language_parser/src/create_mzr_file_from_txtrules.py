#!/usr/bin/python

###############################################################################
# BNGMZRConverter - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
#
# Moleculizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Moleculizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Moleculizer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#                     
#   
###############################################################################

import getopt, sys, re, pdb
import moleculizer


CURRENT_VERSION = 1.0
CURRENT_VERSION_STR = str( CURRENT_VERSION )

def main():
    # Handle command line options
    options = {}
    options["inputfile"] = ""
    options["outputfile"] = "moleculizer-rules-output.mzr"
    options["verbose"] = True

    processOptions(options)

    if not ensureAppropriateOptions(options):
        usage()
        sys.exit(0)

    rulesFile = open(options["inputfile"]).readlines()

    parameterBlock = []
    modificationsBlock  =  []
    molsBlock = []
    allostericPlexes = []
    allostericOmnis = []
    reactionRulesBlock = []
    dimerizationGenBlock = []
    omniGenBlock = []
    uniMolGenBlock = []
    explicitSpeciesBlock = []
    speciesStreamBlock = []

    try:
        parameterBlock, modificationsBlock, molsBlock, allostericPlexes, allostericOmnis, \
        reactionRulesBlock, dimerizationGenBlock, omniGenBlock, uniMolGenBlock, \
        explicitSpeciesBlock, speciesStreamBlock = parseBlockTypesFromRulesFile( rulesFile )
        
    except Exception, e:
        print e.message
        barf("Error, Text rules file could not be parsed.")


    if options["verbose"]:
        print "Parameters:\n\t", parameterBlock
        print "Modifications:\n\t", modificationsBlock
        print "Mols:\n\t", molsBlock
        print "Allosteric Plexes:\n\t", allostericPlexes
        print "allosteric omnis:\n\t", allostericOmnis
        print "Reaction Rules:\n\t", reactionRulesBlock
        print "Dimerization Gens:\n\t", dimerizationGenBlock
        print "OmniGens:\n\t", omniGenBlock
        print "UniMols:\n\t", uniMolGenBlock
        print "ExplicitSpecies:\n\t", explicitSpeciesBlock
        print "Species Streams:\n\t", speciesStreamBlock

    outputFile = moleculizer.MoleculizerRulesFile( options["outputfile"] )

    outputFile.addParameterBlock( parameterBlock )
    outputFile.addModicationsBlock( modificationsBlock )
    outputFile.addMolsBlock( molsBlock )
    outputFile.addAllostericPlexesBlock( allostericPlexes )
    outputFile.addAllostericOmnisBlock( allostericOmnis )
    outputFile.addReactionRulesBlock( reactionRulesBlock, dimerizationGenBlock, \
                                          omniGenBlock, uniMolGenBlock )
    outputFile.addExplicitSpeciesBlock( explicitSpeciesBlock )
    outputFile.addSpeciesStreamsBlock( speciesStreamBlock )

#    outputFile.initialize_DEBUG()

    outputFile.write()
    outputFile.close()
    
    print("Done.")

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
	if opt in ("--version", "-v"):
	    version()
            sys.exit(0)

	if opt in ("--file", "-f"):
	    optionsDict["inputfile"] = atr
	    fileRequired = True
	if opt in ("--output", "-o"):
	    optionsDict["outputfile"] = atr
        if opt in ("--mixin", "-m"):
            optionsDict["events_and_streams_mixin"] = atr

def ensureAppropriateOptions(options):
    if not options["inputfile"]:
        return False
    else:
        return True

def usage():
    print "create_mzr_file_from_txtrules.py version %s" % CURRENT_VERSION
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

    print "Moleculizer is free software; you can redistribute it and/or modify"
    print "it under the terms of the GNU Lesser General Public License as published by"
    print "the Free Software Foundation; either version 3 of the License, or"
    print "(at your option) any later version."
    print
    print "Moleculizer is distributed in the hope that it will be useful,"
    print "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    print "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    print "GNU Lesser General Public License for more details."
    print
    print "You should have received a copy of the GNU Lesser General Public License"
    print "along with Moleculizer; if not, write to the Free Software"
    print "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA"
    print 
    print "Original Author:"
    print "Nathan Addy, Scientific Programmer"
    print "The Molecular Sciences Institute"
    print "                  Email: addy@molsci.org"
    print ""

def parseBlockTypesFromRulesFile(textRulesFile):
    textRulesFile = [re.sub("#.*$", "", x) for x in textRulesFile] # Delete all comments
    # textRulesFile = [re.sub("//.*$", "", x) for x in textRulesFile] # Delete all comments
    textRulesFile = [re.sub(r"\s*", "", x) for x in textRulesFile] # Delete all whitespace
    textRulesFile = [x.strip() for x in textRulesFile] # Strip it for good measure
    textRulesFile = [x for x in textRulesFile if x != ""] # This must be last, because line.strip() results in some empty lines.

# This will attach lines connected with "\"
#     bngFile = '\n'.join(bngFile)
#     bngFile = re.sub(r"\\\s*\n\s*", " ", bngFile)
#     bngFile = bngFile.split("\n")

    parameterBlock = []
    modificationsBlock = []
    molsBlock = []
    allostericPlexes = []
    allostericOmnis = []
    reactionRulesBlock = []
    dimerizationGenBlock = []
    omniGenBlock = [] 
    uniMolGenBlock = [] 
    explicitSpeciesBlock = [] 
    speciesStreamBlock = []

#     textRulesFile = '\n'.join(textRulesFile)
#     textRulesFile = re.sub(r"\\\s*\n\s*", " ", textRulesFile)
#     textRulesFile = textRulesFile.split("\n")

    blockCodes = ["Parameters", "Modifications", "Mols", "Allosteric-Plexes", "Allosteric-Omnis", 
                  "Reaction-Rules", "Dimerization-Gens", "Omni-Gens", "Uni-Mol-Gens", 
                  "Explicit-Species", "Species-Streams" ] 

    blockObjNdx = -1
    blockDataObj = [ (blockCodes[0], parameterBlock), \
                     (blockCodes[1], modificationsBlock), \
                     (blockCodes[2], molsBlock), \
                     (blockCodes[3], allostericPlexes), 
                     (blockCodes[4], allostericOmnis), 
                     (blockCodes[5], reactionRulesBlock), \
                     (blockCodes[6], dimerizationGenBlock), \
                     (blockCodes[7], omniGenBlock), \
                     (blockCodes[8], uniMolGenBlock), \
                     (blockCodes[9], explicitSpeciesBlock),\
                     (blockCodes[10], speciesStreamBlock) ]

    currentDmp = []

    assert( textRulesFile[0].startswith("==="))

    blockObjNdx = -1
    for line in textRulesFile:
        if line.startswith("==="):
            blockObjNdx = returnNewIndex(line, blockObjNdx, blockDataObj)
            currentDmp = blockDataObj[blockObjNdx][1]
        else:
            currentDmp.append(line)


    return getFormattedArray(parameterBlock), getFormattedArray(modificationsBlock), getFormattedArray(molsBlock), getFormattedArray(allostericPlexes), getFormattedArray(allostericOmnis), \
        getFormattedArray(reactionRulesBlock), getFormattedArray(dimerizationGenBlock), getFormattedArray(omniGenBlock), getFormattedArray(uniMolGenBlock), \
        getFormattedArray(explicitSpeciesBlock), getFormattedArray(speciesStreamBlock)

def returnNewIndex(lineOfText, currentNdx, blockObjData):
    key = lineOfText.strip().strip("=").strip()

    while currentNdx < len(blockObjData):
        currentNdx += 1

        if key == blockObjData[currentNdx][0]:
            return currentNdx

    return -1

def barf(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.write("Crashing....\n")
    sys.exit(1)

def printerror(msg):
    sys.stderr.write(msg + '\n')
    return

def getFormattedArray( arrayToFormat ):
    tmpArray = getBalancedArray( arrayToFormat )
    tmpString = "".join( tmpArray )
    if tmpString == "": 
        return []

    assert( tmpString[-1] == ";" )
    tmpArray = tmpString.split(";") 
    tmpArray.pop() # Last entry is blank
    tmpArray = [tok + ";" for tok in tmpArray]
    
    return tmpArray
    

def getBalancedArray( arrayToBalance ):
    if not EachEntryIsParenBalanced( arrayToBalance ):
    # Combine the ..., ndx_i, ndx_(i+1) where ndx_i is the smallest i not balanced
        return getBalancedArray( GetIncrementallyBetterArray( arrayToBalance ) )
    else:
        return arrayToBalance

def GetIncrementallyBetterArray( anArray ):

    values = [ StringIsParenBalenced(x) for x in anArray]

    # This is correct: this function should only be used if the array does not pass
    # EachEntryIsParenBalanced.
    assert( False in values)

    badNdx = values.index( False )
    
    combinedTokens = anArray[badNdx] + anArray[badNdx + 1]

    returnArray = anArray[ : badNdx]
    returnArray.append( combinedTokens )
    returnArray.extend( anArray[badNdx + 2 : ] )
    
    return returnArray
    
def EachEntryIsParenBalanced( array ):
    entries = [ StringIsParenBalenced(x) for x in array ]
    
    returnVal = True
    for val in entries:
        returnVal = returnVal and val

    return returnVal

def StringIsParenBalenced(line):
    return ( line.count("(") == line.count(")") and 
             line.count("[") == line.count("]") and 
             line.count("{") == line.count("}") )

if __name__ == "__main__":
    main()



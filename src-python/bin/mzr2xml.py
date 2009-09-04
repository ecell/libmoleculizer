##!/usr/bin/env/python

###############################################################################
# mzr2xml - Converts rule files in mzr format to intermediate form xml files. 
# Copyright (C) 2007, 2008, 2009 Nathan Addy
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
#   Nathan Addy                  	Email: nathan.addy@gmail.com
#   
###############################################################################

import pdb

import getopt, sys, re
from libmoleculizer import MzrLanguageParser

PROGRAM_NAME = "mzr2xml"
CURRENT_VERSION = 1.1
CURRENT_VERSION_STR = str( CURRENT_VERSION )

def main():
    # Handle command line options
    options = {}
    options["input_file"] = ""
    options["output_file"] = "moleculizer-rules-file.xml"
    options["verbose"] = False

    processOptions(options)

    # There is a mandatory -f or --file= command corresponding to the "input_file"
    if not ensureAppropriateOptions(options):
        usage()
        sys.exit(0)

    rules_file = open(options["input_file"]).readlines()

    parameter_block = []
    modifications_block  =  []
    mols_block = []
    allostericPlexes = []
    allostericOmnis = []
    reactionRules_block = []
    dimerizationGen_block = []
    omniGen_block = []
    explicitSpecies_block = []
    speciesStream_block = []

    try:
        raw_mzr_file = MzrLanguageParser.MzrRawRulesFile(rules_file)

        parameter_block = raw_mzr_file.parameter_block
        modifications_block = raw_mzr_file.modifications_block
        allosteric_plexes = raw_mzr_file.allosteric_plexes_block
        allosteric_omnis = raw_mzr_file.allosteric_omnis_block
        dimerization_gens = raw_mzr_file.association_reactions_block
        omni_gen_block = raw_mzr_file.transformation_reaction_block
        explicit_species = raw_mzr_file.explicit_species_block
        species_streams =  = raw_mzr_file.species_class_block
        
        del(rules_file)
        del(raw_mzr_file)
        
    except Exception, e:
        print e.message
        sys.stderr.write("Error, Text rules file could not be parsed.\n")
        raise e

    if options["verbose"]:
        print("Parameters:\n\t", parameter_block )
        print("Modifications:\n\t", modifications_block)
        print("Molecules:\n\t", mols_block)
        print("Explicit-Allostery:\n\t", allostericPlexes)
        print("Allostery-Classes omnis:\n\t", allostericOmnis)
        print("Association-Reactions:\n\t", dimerizationGen_block)
        print("Transformation-Reactions:\n\t", omniGen_block)
        print("Explicit-Species:\n\t", explicitSpecies_block)
        print("Species-Classes:\n\t", speciesStream_block)

    output_file = moleculizer.MoleculizerRulesFile( options["output_file"] )

    output_file.addParameterBlock( parameterBlock )
    output_file.addModicationsBlock( modificationsBlock )
    output_file.addMolsBlock( molsBlock )
    output_file.addAllostericPlexesBlock( allostericPlexes )
    output_file.addAllostericOmnisBlock( allostericOmnis )
    output_file.addReactionRulesBlock( reactionRulesBlock, dimerizationGenBlock, \
                                          omniGenBlock, [] )
    output_file.addExplicitSpeciesBlock( explicitSpeciesBlock )
    output_file.addSpeciesStreamsBlock( speciesStreamBlock )

#    output_file.initialize_DEBUG()

    output_file.write()
    output_file.close()
    
    print("Done.")

def processOptions( optionsDict ):
    try:
	options, arguments = getopt.getopt(sys.argv[1:], "hf:o:vm:", ["help", "file=", "output=", "verbose", "version", "mixin="])
    except getopt.error, msg:
	print msg
	print "For help, use --help"
	sys.exit(0)

    for opt, atr in options:

	if opt in ("--help", "-h"):
	    usage()
            sys.exit(0)
	if opt in ("--version",):
	    version()
            sys.exit(0)

        if opt in ("--verbose", "-v"):
            optionsDict["verbose"] = True
	if opt in ("--file", "-f"):
	    optionsDict["input_file"] = atr
	    fileRequired = True
	if opt in ("--output", "-o"):
	    optionsDict["output_file"] = atr
        if opt in ("--mixin", "-m"):
            optionsDict["events_and_streams_mixin"] = atr

def ensureAppropriateOptions(options):
    if not options["input_file"]:
        return False
    else:
        return True

def usage():
    print "%s version %s" % (PROGRAM_NAME, CURRENT_VERSION)
    print "\t-- Creates a Moleculizer XML input file from a human-readable text version of biochemical network rules."
    print
    print "Usage:"
    print "mzrrulestoxmlconverter.py -f inputfile [-o outputfile] [-m eventsandstreamsmixin]"
    print "mzrrulestoxmlconverter.py [-h] or [--help]"
    print 
    print "Author:"
    print "Nathan Addy (nathan.addy@gmail.com)"
    print "Molecular Sciences Institute, Berkeley, CA"
    print "March 3, 2009"
    print
    print 
    print "This is a parser that takes in a moleculizer rules file in text file format and writes it to a moleculizer xml file"
    print "for input into a moleculizer run."

def version():
    print "%s %s -- Copyright (C) 2009 The Moleculizer Sciences Institute" % (PROGRAM_NAME, CURRENT_VERSION)
    print
    print "Author:"
    print "Nathan Addy (nathan.addy@gmail.com)"
    print "Molecular Sciences Institute, Berkeley, CA"
    print "March 3, 2009"
    print 
    print "Licensed under LGPL v3."
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
    explicitSpeciesBlock = [] 
    speciesStreamBlock = []

#     textRulesFile = '\n'.join(textRulesFile)
#     textRulesFile = re.sub(r"\\\s*\n\s*", " ", textRulesFile)
#     textRulesFile = textRulesFile.split("\n")

    blockCodes = ["Parameters", "Modifications", "Molecules", "Explicit-Allostery", "Allostery-Classes", 
                  "Reaction-Rules", "Association-Reactions", "Transformation-Reactions", 
                  "Explicit-Species", "Species-Classes" ] 

    blockObjNdx = -1
    blockDataObj = [ (blockCodes[0], parameterBlock), \
                     (blockCodes[1], modificationsBlock), \
                     (blockCodes[2], molsBlock), \
                     (blockCodes[3], allostericPlexes), 
                     (blockCodes[4], allostericOmnis), 
                     (blockCodes[5], reactionRulesBlock), \
                     (blockCodes[6], dimerizationGenBlock), \
                     (blockCodes[7], omniGenBlock), \
                     (blockCodes[8], explicitSpeciesBlock),\
                     (blockCodes[9], speciesStreamBlock) ]

    currentDmp = []

    assert( textRulesFile[0].startswith("="))

    blockObjNdx = -1
    for line in textRulesFile:
        if line.startswith("="):
            blockObjNdx = returnNewIndex(line, blockDataObj)
            currentDmp = blockDataObj[blockObjNdx][1]
        else:
            currentDmp.append(line)


    return getFormattedArray(parameterBlock), getFormattedArray(modificationsBlock), getFormattedArray(molsBlock), getFormattedArray(allostericPlexes), getFormattedArray(allostericOmnis), \
        getFormattedArray(reactionRulesBlock), getFormattedArray(dimerizationGenBlock), getFormattedArray(omniGenBlock), \
        getFormattedArray(explicitSpeciesBlock), getFormattedArray(speciesStreamBlock)

def returnNewIndex(lineOfText, blockObjData):
    key = lineOfText.strip().strip("=").strip()

    for ndx in range(len(blockObjData)):
        if key == blockObjData[ndx][0]:
            return ndx
    raise Exception("Section title '%s' cannot be found" % key)

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

    try:
        assert( tmpString[-1] == ";" )
    except:
        raise Exception("Error parsing block '%s'.  Line does not end in ';'." % repr(arrayToFormat))

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
    pdb.set_trace()
    main()



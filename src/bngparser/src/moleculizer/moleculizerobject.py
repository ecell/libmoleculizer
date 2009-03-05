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

import re, pdb
import util

import mzrexceptions as MzrExceptions

from expressionevaluation import SymbolicExpressionEvaluator
from moleculeinterpreter import MoleculeDictionary
from seedspeciesinterpreter import SeedSpecies
from reactionrulesinterpreter import ReactionRules
from xmlobject import XmlObject
from moleculizermol import MoleculizerMol

class MoleculizerObject:
    def __init__(self, filename):
        
        # Self explanatory...
        self.outputFileName = filename

        # The data which will be parsed and used to create the 
        # moleculizer file.
        self.parameterBlock = []
        self.molsBlock = []
        self.allostericPlexes = []
        self.allostericOmnis = []
        self.reactionRulesBlock = []
        self.dimerizationGenBlock = []
        self.omniGenBlock = []
        self.uniMolGenBlock = []
        self.explicitSpeciesBlock = []
        self.speciesStreamBlock = []

        # These are the objects that will be used to process the parsed 
        # data.

	self.parameterEE = 0 
	self.moleculeDictionary = 0
        self.allostericData = 0
	self.reactionRules = 0
        self.explicitSpecies = 0
        self.speciesStreams = 0

    def getOutputFileName(self):
        return self.outputFileName

    def write(self):
    	# This creates objects that parse and analyse the reactionBlocks.
	self.__processData()

	# This lets each of the objects write its xml output to the master xml object.
        self.openXmlFile = open(self.outputFileName, 'w')
	self.__writeOutput(self.openXmlFile)
	return

    def initialize_DEBUG(self):
        self.__processData()

    def close(self):
        self.openXmlFile.close()

    def addParameterBlock(self, parameterBlock, overwrite = False):
	if self.parameterBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a parameter block twice.")
	self.parameterBlock = parameterBlock[:]

    def addMolsBlock(self, molsBlock):
        if self.molsBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a mols block twice.")
        self.molsBlock = molsBlock[:]

    def addAllostericPlexesBlock(self, apBlock, overwrite = False):
        if self.allostericPlexes and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add an allosteric plexes block twice.")
        self.allostericPlexes = apBlock[:]

    def addAllostericOmnisBlock(self, aoBlock, overwrite = False):
        if self.allostericOmnis and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add an allosteric omnis block twice.")
        self.allostericOmnis = aoBlock[:]

    def addReactionRulesBlock( self, rrBlock, overwrite = False):
        if self.reactionRulesBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a reaction rules block twice.")
        self.reactionRulesBlock = rrBlock[:]

    def addDimerizationGensBlock( self, dgBlock, overwrite = False):
        if self.dimerizationGenBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a dimerization gen block twice.")
        self.dimerizationGenBlock = dgBlock[:]

    def addOmniGensBlock(self, ogBlock, overwrite = False):
        if self.omniGenBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a omni gen block twice.")
        self.omniGenBlock=ogBlock[:]

    def addUniMolGensBlock( self, umBlock, overwrite = False):
        if self.uniMolGenBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a unimol block twice.")
        self.uniMolGenBlock = umBlock[:]

    def addExplicitSpeciesBlock( self, esBlock, overwrite = False):
        if self.explicitSpeciesBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add an explicit species block twice.")
        self.explicitSpeciesBlock = esBlock[:]

    def addSpeciesStreamsBlock(self, ssBlock, overwrite = False):
        if self.speciesStreamBlock and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add a species stream block twice.")
        self.speciesStreamBlock = ssBlock[:]

    def __processData(self):
        ## "And for what you are about to see next, we must enter, quietly, into the
        ##  realm of genius..."

	self.parameterEE = self.__processParameterBlock(self.parameterBlock)
 	self.moleculeDictionary = self.__processMolsBlock(self.molsBlock )
        self.allostericData = self.__processAllostericRulesBlocks(self.allostericPlexes, self.allostericOmnis)
        self.reactionRules = self.__processReactionRulesBlocks( self.reactionRulesBlock, self.dimerizationGenBlock, self.omniGenBlock, self.uniMolGenBlock)
        self.explicitSpecies = self.__processExplicitSpeciesBlock( self.explicitSpeciesBlock)
        self.speciesStreams = self.__processSpeciesStreamBlock( self.speciesStreamBlock)

        return


    def __processParameterBlock(self, parameterBlock):
        ## This function works by instantiating a SymbolicExpressionEvaluator with a
        ## parameterBlock and returns the object, which can then be passed to other objects.
        ## It's "evaluateExpression" function will take an expression and return a numerical
        ## value.

	tmpSEE = SymbolicExpressionEvaluator(self.parameterBlock)
	return tmpSEE

    def __processMolsBlock(self, moleculeBlock):
        ## This function works by taking instantiating a MoleculeDictionary with a
        ## moleculeDescriptionBlock and a SymbolicExpressionEvaluator and returns the object,
        ## which can then be passed to other objects.
        assert(self.parameterEE)

	moleculeDictionary = MoleculeDictionary(moleculeBlock, self.parameterEE)
        return moleculeDictionary

    def __processAllostericRulesBlocks( self, allostericPlexBlock, allostericOmniBlock):
        return 0

    def __processReactionRulesBlocks( self, rxnRulesBlock, dimerBlock, omniGenBlock, uniGenBlock):
        return 0

    def __processExplicitSpeciesBlock( self, explicitSpeciesBlock):
        return 0

    def __processSpeciesStreamBlock( self, ssBlock):
        return 0

    
#     def __processMoleculeTypeBlock(self, moleculeBlock):
#         ## This function works by taking instantiating a MoleculeDictionary with a
#         ## moleculeDescriptionBlock and a SymbolicExpressionEvaluator and returns the object,
#         ## which can then be passed to other objects.
#         assert(self.parameterEE)

# 	moleculeDictionary = MoleculeDictionary(moleculeBlock, self.parameterEE)
#         return moleculeDictionary

#     def __processSeedSpeciesBlock(self, seedSpeciesBlock):
#         ## This works by instantiating a SeedSpecies object, by passing it a seed species block,
#         ## a MoleculeDictionary, and a SymbolicExpressionEvaluator.  It will return a
#         ## seedSpecies object.

#         assert(self.parameterEE)
#         assert(self.moleculeDictionary)
# 	seedSpeciesObj = SeedSpecies(seedSpeciesBlock, self.moleculeDictionary, self.parameterEE)

#         return seedSpeciesObj

#     def __processReactionRulesBlock(self, rxnRulesBlock):
#         ## This does the same thing: pass it a reaction rules block, a MoleculeDictionary,
#         ## and a SymbolicExpressionEvaluator and it will return a ReactionRules object.

#         assert(self.parameterEE)
#         assert(self.moleculeDictionary)
        
# 	reactionRules = ReactionRules(rxnRulesBlock, self.moleculeDictionary, self.parameterEE)
#         return reactionRules








    def __writeOutput(self, openXMLFile):
        xmlobject = self.__constructXMLRepresentation()
        xmlobject.writeall(openXMLFile)

    def __constructXMLRepresentation(self):
        rootNode = XmlObject("moleculizer-input")
        modelElmt = XmlObject("model")
        modelElmt.attachToParent(rootNode)

        self.__addModifications( modelElmt )
        self.__addMols( modelElmt )
        self.__addAllostericPlexes( modelElmt )
        self.__addAllostericOmnis( modelElmt )
        self.__addReactionGens( modelElmt )
        self.__addExplicitSpecies( modelElmt )
        self.__addExplicitReactions( modelElmt )
        self.__addVolume( modelElmt )

        # This adds the streams and events elements.
        if self.streamsAndEventsFile:
            rootNode.addSubFile(self.streamsAndEventsFile)

        return rootNode

    def __addModifications(self, parentObject):
        modificationListElmt = XmlObject("modifications")
        modificationListElmt.attachToParent(parentObject)
        self.moleculeDictionary.addModifications(modificationListElmt)

    def __addMols(self, parentObject):
        molsElmt = XmlObject("mols")
        molsElmt.attachToParent(parentObject)
        self.moleculeDictionary.addMols(molsElmt)

    def __addAllostericPlexes(self, parentObject):
        allostericPlexesElmt = XmlObject("allosteric-plexes")
        allostericPlexesElmt.attachToParent(parentObject)

    def __addAllostericOmnis(self, parentObject):
        allostericOmnisElmt = XmlObject("allosteric-omnis")
        allostericOmnisElmt.attachToParent(parentObject)
        self.reactionRules.addAllostericOmnis(allostericOmnisElmt)

    def __addReactionGens(self, parentObject):
        rxnGenElmt = XmlObject("reaction-gens")
        rxnGenElmt.attachToParent(parentObject)
        self.reactionRules.writeReactionGensElement(rxnGenElmt)

    def __addExplicitSpecies(self, parentObject):
        explicitSpeciesElmt = XmlObject("explicit-species")
        explicitSpeciesElmt.attachToParent(parentObject)
        self.seedSpecies.writeExplicitSpecies(explicitSpeciesElmt)

    def __addExplicitReactions(self, parentObject):
        explicitRxnElmt = XmlObject("explicit-reactions")
        parentObject.addSubElement(explicitRxnElmt)

    def __addVolume(self, parentObject):
        ## So far, there is a parameter "Cell_volume" that does this.
        volumeElement = XmlObject("volume")
        volumeElement.addAttribute("liters", self.parameterEE.getVariableValue("Cell_volume"))
        parentObject.addSubElement(volumeElement)


    
	


    


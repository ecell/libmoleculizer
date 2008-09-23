###############################################################################
# BNGMZRConverter - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007, 2008 The Molecular Sciences Institute
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
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Contact information:
#   Nathan Addy, Scientific Programmer	Voice: 510-981-8748
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#   2168 Shattuck Ave.                  
#   Berkeley, CA 94704
###############################################################################

import re, pdb
import util

import mzrexceptions as MzrExceptions

from expressionevaluation import AlphaWikiSymbolicExpressionEvaluator
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
	self.moleculeTypeBlock = []
	self.seedSpeciesBlock = []
	self.reactionRuleBlock = []
        self.streamsAndEventsFile = 0

        # These are the objects that will be used to process the parsed 
        # data.
	self.parameterEE = 0 
	self.moleculeDictionary = 0
	self.seedSpecies = 0
	self.reactionRules = 0

    def getOutputFileName(self):
        return self.outputFileName

    def write(self):
    	# This creates objects that parse and analyse the reactionBlocks.
	self.__processData()

	# This lets each of the objects write its xml output to the master xml object.
        self.openXmlFile = open(self.outputFileName, 'w')
	self.__writeOutput(self.openXmlFile)
	return

    def close(self):
        self.openXmlFile.close()

    def addParameterBlock(self, parameterBlock):
	if self.parameterBlock: raise MzrExceptions.MoleculizerException("Error: Cannot add a parameter block twice.")
	self.parameterBlock = parameterBlock[:]

    def addMoleculeTypeBlock(self, moleculeTypeBlock):
	if self.moleculeTypeBlock: raise MzrExceptions.MoleculizerException("Error: Cannot add a molecule types block twice")
	self.moleculeTypeBlock = moleculeTypeBlock[:]

    def addSeedSpeciesBlock(self, seedSpeciesBlock):
	if self.seedSpeciesBlock: raise MzrExceptions.MoleculizerException("Error:  Cannot add a Seed Species block twice.")
	self.seedSpeciesBlock = seedSpeciesBlock[:]

    def addReactionRuleBlock(self, reactionRuleBlock):
	if self.reactionRuleBlock: raise MzrExceptions.MoleculizerException("Error: Cannot add a Reaction Rule Block twice.")
	self.reactionRuleBlock = reactionRuleBlock[:]

    def addStreamsAndEventsMixin(self, streamsAndEventsFile):
        self.streamsAndEventsFile = streamsAndEventsFile

    def __processData(self):
        ## "And for what you are about to see next, we must enter, quietly, into the
        ##  realm of genius..."

	self.parameterEE = self.__processParameterBlock(self.parameterBlock)
 	self.moleculeDictionary = self.__processMoleculeTypeBlock(self.moleculeTypeBlock)
        self.seedSpecies = self.__processSeedSpeciesBlock(self.seedSpeciesBlock)
        self.reactionRules = self.__processReactionRulesBlock(self.reactionRuleBlock)

        return

    def __processParameterBlock(self, parameterBlock):
        ## This function works by instantiating a SymbolicExpressionEvaluator with a
        ## parameterBlock and returns the object, which can then be passed to other objects.
        ## It's "evaluateExpression" function will take an expression and return a numerical
        ## value.

	tmpSEE = AlphaWikiSymbolicExpressionEvaluator(self.parameterBlock)
	return tmpSEE
    
    def __processMoleculeTypeBlock(self, moleculeBlock):
        ## This function works by taking instantiating a MoleculeDictionary with a
        ## moleculeDescriptionBlock and a SymbolicExpressionEvaluator and returns the object,
        ## which can then be passed to other objects.
        assert(self.parameterEE)

	moleculeDictionary = MoleculeDictionary(moleculeBlock, self.parameterEE)
        return moleculeDictionary

    def __processSeedSpeciesBlock(self, seedSpeciesBlock):
        ## This works by instantiating a SeedSpecies object, by passing it a seed species block,
        ## a MoleculeDictionary, and a SymbolicExpressionEvaluator.  It will return a
        ## seedSpecies object.

        assert(self.parameterEE)
        assert(self.moleculeDictionary)
	seedSpeciesObj = SeedSpecies(seedSpeciesBlock, self.moleculeDictionary, self.parameterEE)

        return seedSpeciesObj

    def __processReactionRulesBlock(self, rxnRulesBlock):
        ## This does the same thing: pass it a reaction rules block, a MoleculeDictionary,
        ## and a SymbolicExpressionEvaluator and it will return a ReactionRules object.

        assert(self.parameterEE)
        assert(self.moleculeDictionary)
        
	reactionRules = ReactionRules(rxnRulesBlock, self.moleculeDictionary, self.parameterEE)
        return reactionRules

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


    
	


    


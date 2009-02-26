###############################################################################
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    
#
###############################################################################

import pdb

import re

import util
from moleculizer_xcpt import *
from xmlobject import XmlObject

class MoleculizerRulesFile:
    """
    This object acts as an parsing thing that outputs moleculizer files xml,
    suitable for processing internally by a mzr::moleculizer instance."""

    @staticmethod
    def block_posses_sanity_checks( linearray ):
        everyLineEndsWithSemiColon = [ x[-1] == ";" for x in linearray]
        noWhiteSpace = [ (x.count("\n") + x.count(" ") + x.count("\t") == 0) for x in linearray]

        return reduce(util.And, everyLineEndsWithSemiColon) and reduce(util.And, noWhiteSpace)

    @staticmethod
    def snip_block_lines( linearray ):
        everyLineEndsWithSemiColon = [ x[-1] == ";" for x in linearray]
        assert( reduce(util.And, everyLineEndsWithSemiColon) )

        return [x[:-1] for x in linearray]
    

    def __init__(self, filename):

        self.outputFileName = filename

        # These are the lines of input, in one statement per line form 
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

        # A section is an intermediate between a rules file (they have lines, for example,
        # and can answer questions about what has been parsed ) and an xml section (it can
        # write out an xml section -
        # Parameters doesn't write anything out currently, but easily could
        # moleculeSection --- <modifications> and <mols>
        # allosterySection --- <allosteric-plexes> and <allosteric-omnis>
        # reactionRulesSection -- <reaaction-gens>?
        # explicitSpeciesSection -- <explicit-species>
        # speciesStreamSection -- <species-streams>
	self.parameterSection = 0 
	self.moleculeSection = 0
        self.allosterySection = 0
	self.reactionRulesSection = 0
        self.explicitSpeciesSection = 0
        self.speciesStreamSection = 0

    # output_file = property( getOutputFileName )
    @property
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

    def initialize_DEBUG(self):
        self.__processBlocksData()

    def addParameterBlock(self, parameterBlock, overwrite = False):
	if self.parameterBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a parameter block twice.")

        if not self.block_passes_sanity_checks( parameterBlock ):
            raise InsaneBlockOnTheLooseException(parameterBlock, "parameter block")

	self.parameterBlock = parameterBlock[:]
        self.parameterEE = SymbolicExpressionEvaluator( self.parameterBlock )

    def addMolsBlock(self, molsBlock):
        if self.molsBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a mols block twice.")

        if not self.block_passes_sanity_checks( molsBlock ):
            raise InsaneBlockOnTheLooseException(molsBlock, "mols block")

        self.molsBlock = molsBlock[:]
        self.molsSection = MolsSection( moleculeBlock )

    def addAllostericPlexesBlock(self, apBlock, overwrite = False):
        if self.allostericPlexes and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add an allosteric plexes block twice.")

        if not self.block_passes_sanity_checks( apBlock ):
            raise InsaneBlockOnTheLooseException(apBlock, "allosteric plexes block")

        self.allostericPlexes = apBlock[:]
        self.allostericPlexesSection = AllostericPlexes( self.allostericPlexes )

    def addAllostericOmnisBlock(self, aoBlock, overwrite = False):
        if self.allostericOmnis and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add an allosteric omnis block twice.")
        
        if not self.block_passes_sanity_checks( aoBlock ):
            raise InsaneBlockOnTheLooseException( aoBlock, "allosteric omnis block")
        self.allostericOmnis = aoBlock[:]
        self.allostericOmnisSection = AllostericOmnis( self.allostericOmnis )

    def addReactionRulesBlock( self, rrBlock, overwrite = False):
        if self.reactionRulesBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a reaction rules block twice.")


        if not self.block_passes_sanity_checks( rrBlock ):
            raise InsaneBlockOnTheLooseException(rrBlock, "reaction rules")

        if not self.block_passes_sanity_checks( dgBlock ):
            raise InsaneBlockOnTheLooseException(dgBlock, "dimerization gen block")

        if not self.block_passes_sanity_checks( omniGenBlock ):
            raise InsaneBlockOnTheLooseException(omniGenBlock, "omni-gen block")

        if not self.block_passes_sanity_checks( uniMolGenBlock ):
            raise InsaneBlockOnTheLooseException(uniMolGenBlock, "uni-mol-gen block")

        self.reactionRulesBlock = rrBlock[:]
        self.dimerizationGenBlock = dgBlock[:]
        self.omniGenBlock=ogBlock[:]
        self.uniMolGenBlock = umBlock[:]


        self.reactionGensSection = ReactionGens( self.reactionRulesBlock, 
                                                 self.dimerizationGenBlock, 
                                                 self.omniGenBlock, 
                                                 self.uniMolGenBlock,
                                                 self.molsSections)

    def addExplicitSpeciesBlock( self, esBlock, overwrite = False):
        if self.explicitSpeciesBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add an explicit species block twice.")
        
        if not self.block_passes_sanity_checks( esBlock  ):
            raise InsaneBlockOnTheLooseException(esBlock, "explicit-species")
        
        self.explicitSpeciesBlock = esBlock[:]
        self.explicitSpeciesSection = ExplicitSpecies( self.explicitSpecies )


    def addSpeciesStreamsBlock(self, ssBlock, overwrite = False):
        if self.speciesStreamBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a species stream block twice.")
        
        if not self.block_passes_sanity_checks( ssBlock ):
            raise InsaneBlockOnTheLooseException(ssBlock, "")
        
        self.speciesStreamBlock = ssBlock[:]
        self.speciesStreamSection = SpeciesStream( self.speciesStreamBlock )

    def __processBlockData(self):
        # The parameterEE (parameter expression evaluator can report on calculated values
        # from the parameters section, as well as evaluate mathematical expressions using 
        # lazy evaluation.

	self.parameterEE = SymbolicExpressionEvaluator( self.parameterBlock )
        MoleculizerSection.setParameterEvaluator( self.parameterEE )

 	self.molsSection = MolsSection( moleculeBlock )

        self.reactionGensSection = ReactionGens( self.reactionRulesBlock, self.dimerizationGenBlock, \
                                                 self.omniGenBlock, self.uniMolGenBlock, \
                                                 self.molsSections)

        self.explicitSpeciesSection = ExplicitSpecies( self.explicitSpecies, self.mols )
        
        self.speciesStreams = SpeciesStreams( self.speciesStreams, self.mols)
        return


#     def __processAllostericRulesBlocks( self, allostericPlexBlock, allostericOmniBlock):
#         return 0

#     def __processReactionRulesBlocks( self, rxnRulesBlock, dimerBlock, omniGenBlock, uniGenBlock):
#         return 0

#     def __processExplicitSpeciesBlock( self, explicitSpeciesBlock):
#         return 0

#     def __processSpeciesStreamBlock( self, ssBlock):
#         return 0
    
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


    
	


    


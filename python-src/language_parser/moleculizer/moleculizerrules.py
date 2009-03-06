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

from xmlobject import XmlObject

from sectionparameter import SymbolicExpressionEvaluator
from sectionmodifications import ModificationsSection
from sectionmols import MolsSection
from sectionallostery import AllostericPlexesSection, AllostericOmnisSection
from sectionreactionrules import ReactionRulesSection
from sectionspeciesstreams import SpeciesStreamsSection
from sectionexplicitspeciesblock import ExplicitSpeciesSection

from moleculizer_xcpt import *
class MoleculizerRulesFile:
    """
    This object acts as an parsing thing that outputs moleculizer files xml,
    suitable for processing internally by a mzr::moleculizer instance."""

    @staticmethod
    def BlockPassesSanityCheck( linearray ):
        linearray = [x for x in linearray if x.strip() != ""]
        if len(linearray) == 0: return True
        
        everyLineEndsWithSemiColon = [ x[-1] == ";" and x.count(";") == 1for x in linearray]
        noWhiteSpace = [ (x.count("\n") + x.count(" ") + x.count("\t") == 0) for x in linearray]
        return reduce(util.And, everyLineEndsWithSemiColon) and reduce(util.And, noWhiteSpace)

    def DEBUGPRINT(self):
        print "Hello, from moleculizer!!!!!!"

    def __init__(self, filename):

        self.outputFileName = filename

        # These are the lines of input, in one statement per line form, with no whitespace
        self.parameterBlock = []
        self.modificationsBlock = []
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
	self.parameterSection = 0 
        self.modificationsSection = 0
	self.molsSection = 0
        self.allostericPlexesSection = 0
        self.allostericOmnisSection = 0
	self.reactionGensSection = 0
        self.explicitSpeciesSection = 0
        self.speciesStreamSection = 0

    # output_file = property( getOutputFileName )
    @property
    def getOutputFileName(self):
        return self.outputFileName
        
    def write(self):
        self.openXmlFile = open(self.outputFileName, 'w')

	self.__writeOutput(self.openXmlFile)

	return

    def close(self):
        self.openXmlFile.close()

    def addParameterBlock(self, parameterBlock, overwrite = False):
	if self.parameterBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a parameter block twice.")

        if not self.BlockPassesSanityCheck( parameterBlock ):
            raise InsaneBlockOnTheLooseException(parameterBlock, "parameter block")

	self.parameterBlock = parameterBlock[:]
        self.parameterEE = SymbolicExpressionEvaluator( self.parameterBlock )

    def addModicationsBlock(self, modificationsBlock, overwrite = False):
        if self.modificationsBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a modifications block twice.")

        if not self.BlockPassesSanityCheck( modificationsBlock ):
            raise InsaneBlockOnTheLooseException(modificationsBlock, "modifications block")

        self.modificationsBlock = modificationsBlock[:]
        self.modificationsSection = ModificationsSection( self.modificationsBlock  )

        return


    def addMolsBlock(self, molsBlock):
        if self.molsBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a mols block twice.")

        if not self.BlockPassesSanityCheck( molsBlock ):
            raise InsaneBlockOnTheLooseException(molsBlock, "mols block")

        self.molsBlock = molsBlock[:]
        self.molsSection = MolsSection( molsBlock  )


    def addAllostericPlexesBlock(self, apBlock, overwrite = False):
        if self.allostericPlexes and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add an allosteric plexes block twice.")

        if not self.BlockPassesSanityCheck( apBlock ):
            raise InsaneBlockOnTheLooseException(apBlock, "allosteric plexes block")

        self.allostericPlexes = apBlock[:]
        self.allostericPlexesSection = AllostericPlexesSection( self.allostericPlexes )

    def addAllostericOmnisBlock(self, aoBlock, overwrite = False):
        if self.allostericOmnis and not overwrite: raise MzrExceptions.MoleculizerException("Error: Cannot add an allosteric omnis block twice.")
        
        if not self.BlockPassesSanityCheck( aoBlock ):
            raise InsaneBlockOnTheLooseException( aoBlock, "allosteric omnis block")
        self.allostericOmnis = aoBlock[:]
        self.allostericOmnisSection = AllostericOmnisSection( self.allostericOmnis )

    def addReactionRulesBlock( self, rrBlock, dimerGenBlock, omniGenBlock, uniMolGenBlock, overwrite = False):

        if self.reactionRulesBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a reaction rules block twice.")

        if not self.BlockPassesSanityCheck( rrBlock ):
            raise InsaneBlockOnTheLooseException(rrBlock, "reaction rules")

        if not self.BlockPassesSanityCheck( dimerGenBlock ):
            raise InsaneBlockOnTheLooseException(dimerGenBlock, "dimerization gen block")

        if not self.BlockPassesSanityCheck( omniGenBlock ):
            raise InsaneBlockOnTheLooseException(omniGenBlock, "omni-gen block")

        if not self.BlockPassesSanityCheck( uniMolGenBlock ):
            raise InsaneBlockOnTheLooseException(uniMolGenBlock, "uni-mol-gen block")


        self.reactionRulesBlock.extend( rrBlock )
        self.dimerizationGenBlock.extend( dimerGenBlock )
        self.omniGenBlock.extend( omniGenBlock )
        self.uniMolGenBlock.extend( uniMolGenBlock )

        self.reactionGensSection = ReactionRulesSection( self.reactionRulesBlock, 
                                                         self.dimerizationGenBlock, 
                                                         self.omniGenBlock, 
                                                         self.uniMolGenBlock)

    def addExplicitSpeciesBlock( self, esBlock, overwrite = False):
        if self.explicitSpeciesBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add an explicit species block twice.")
        
        if not self.BlockPassesSanityCheck( esBlock  ):
            raise InsaneBlockOnTheLooseException(esBlock, "explicit-species")
        
        self.explicitSpeciesBlock = esBlock[:]
        self.explicitSpeciesSection = ExplicitSpeciesSection( esBlock )


    def addSpeciesStreamsBlock(self, ssBlock, overwrite = False):
        if self.speciesStreamBlock and not overwrite: 
            raise MzrExceptions.MoleculizerException("Error: Cannot add a species stream block twice.")

        if not self.BlockPassesSanityCheck( ssBlock ):
            raise InsaneBlockOnTheLooseException(ssBlock, "")
        
        self.speciesStreamBlock = ssBlock[:]
        self.speciesStreamSection = SpeciesStreamsSection( ssBlock  )


    def __processAllostericRulesBlocks( self, allostericPlexBlock, allostericOmniBlock):
        return 0
    def __processReactionRulesBlocks( self, rxnRulesBlock, dimerBlock, omniGenBlock, uniGenBlock):
        return 0

    def __processExplicitSpeciesBlock( self, explicitSpeciesBlock):
        return 0

    def __processSpeciesStreamBlock( self, ssBlock):
        return 0
    
    def __writeOutput(self, openXMLFile):
        xmlobject = self.__constructXMLRepresentation()
        xmlobject.writeall(openXMLFile)

    def __constructXMLRepresentation(self):
        
        rootNode = XmlObject("moleculizer-input")
        modelElmt = XmlObject("model")
        modelElmt.attachToParent(rootNode)

        streamsElmt = XmlObject("streams", rootNode)

        self.__addModifications( modelElmt )

        self.__addMols( modelElmt )
        self.__addAllostericPlexes( modelElmt )
        self.__addAllostericOmnis( modelElmt )
        self.__addReactionGens( modelElmt )
        self.__addExplicitSpecies( modelElmt )
        self.__addExplicitReactions( modelElmt )

        self.__addSpeciesStreams( streamsElmt )

        return rootNode

    def __addModifications(self, parentObject):
        # Write me!!!
        modificationsSection = XmlObject("modifications", parentObject)
        
        if self.modificationsSection:
            self.modificationsSection.writeModificationsSections( modificationsSection )
        return

    def __addMols(self, parentObject):
        molsSection = XmlObject("mols", parentObject)

        if self.molsSection:
            self.molsSection.writeMolsSection( molsSection)
        
        return 

    def __addAllostericPlexes(self, parentObject):
        allostericPlexes = XmlObject("allosteric-plexes", parentObject)

        if self.allostericPlexesSection:
            self.allostericPlexesSection.writeAllostericPlexesSection(allostericPlexes)

        return 

    def __addAllostericOmnis(self, parentObject):
        allostericOmnis = XmlObject("allosteric-omnis", parentObject)

        if self.allostericOmnisSection:
            self.allostericOmnisSection.writeAllostericOmnisSection( allostericOmnis )

        return 

    def __addReactionGens(self, parentObject):
        reactionGenElmt = XmlObject("reaction-gens", parentObject)

        
        
        if self.reactionGensSection:
            self.reactionGensSection.writeReactionGensSection( reactionGenElmt )

        return 

    def __addSpeciesStreams( self, parentObject):

        speciesStreamsElement = XmlObject("species-streams", parentObject)

        if self.speciesStreamSection:
            self.speciesStreamSection.writeSpeciesStreamSection( speciesStreamsElement )

    def __addExplicitSpecies(self, parentObject):
        explicitSpeciesElmt = XmlObject("explicit-species", parentObject)
        
        if self.explicitSpeciesSection:
            self.explicitSpeciesSection.writeExplicitSpeciesSection( explicitSpeciesElmt )

        return

    def __addExplicitReactions( self, modelElmt ):
        explicitReactionsElmt = XmlObject("explicit-reactions", modelElmt)
        return 




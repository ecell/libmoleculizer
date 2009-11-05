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

from xmlobject import XmlObject
from bngspeciesdefinition import BngSpeciesInstantiation
import pdb

class SeedSpecies:
    def __init__(self, seedSpeciesBlock, moleculeDictionary, parameterDictionary):
        self.rawSeedSpeciesBlock = seedSpeciesBlock[:]
        self.moleculeDictionary = moleculeDictionary
        self.parameterDictionary = parameterDictionary
        self.seedSpeciesInModel = []

        self.initialize()

    def initialize(self):
        for line in self.rawSeedSpeciesBlock:
            self.parseSeedSpeciesLine(line)

    def parseSeedSpeciesLine(self, seedSpeciesLine):
	seedSpeciesLineComponents = seedSpeciesLine.split()
        assert( len(seedSpeciesLineComponents ) == 2)

        seedSpeciesDef, number = seedSpeciesLineComponents	    
        realNumber = self.parameterDictionary.evaluateExpression(number)

	self.seedSpeciesInModel.append( (BngSpeciesInstantiation(seedSpeciesDef), realNumber) )

    def writeExplicitSpecies(self, parentElmt):

        for ss in self.seedSpeciesInModel:
            plexSpeciesElmt = XmlObject("plex-species")
            plexSpeciesElmt.attachToParent(parentElmt)
            
            speciesDefinition, number = ss

            # Inputs to the species numbers are given in micromolar...
            # The inputs we have been given are in the form: Dig1_tot_conc=Dig1_num*1e6/(Cell_volume*Avogadros_number)
            # So to get num we multiply by (Cell_v * Avog) / 1e6.
            number = number * (6.0221415e23 * 40e-15) / 1e6
            plexSpeciesElmt.addAttribute("name", speciesDefinition.constructEntireComplexInstanceName())

            plexElmt = XmlObject("plex")
            plexElmt.attachToParent(plexSpeciesElmt)

            populationElmt = XmlObject("population")
            populationElmt.attachToParent(plexSpeciesElmt)
            populationElmt.addAttribute("count", int(number))

            # Make this a method of class bndSpeciesDefinition.
            for molNdx in speciesDefinition.moleculeRefsInComplex:
                molReferenceName = speciesDefinition.moleculeRefsInComplex[molNdx]
                molInstanceName = speciesDefinition.moleculeInstancesInComplex[molNdx]
                
                molInstanceElmt = XmlObject("mol-instance")
                molInstanceElmt.addAttribute("name", molInstanceName)
                molReferenceElmt = XmlObject("mol-ref")
                molReferenceElmt.addAttribute("name", molReferenceName)
                molReferenceElmt.attachToParent(molInstanceElmt).attachToParent(plexElmt)

            for bindingNdx in speciesDefinition.bindings:
                bindingElmt = XmlObject("binding")
                bindingElmt.attachToParent(plexElmt)

                for hb in speciesDefinition.bindings[bindingNdx]:
                    molInstanceRefElmt = XmlObject("mol-instance-ref")
                    molInstanceRefElmt.addAttribute("name", speciesDefinition.moleculeInstancesInComplex[hb[0]])
                    bindingSiteRefElmt = XmlObject("binding-site-ref")
                    bindingSiteRefElmt.addAttribute("name", hb[1])
                    bindingSiteRefElmt.attachToParent(molInstanceRefElmt).attachToParent(bindingElmt)

                

                



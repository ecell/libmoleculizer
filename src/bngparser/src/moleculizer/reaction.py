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

import pdb
import util
class Reaction:

    volume = 0.0
    
    def __init__(self, reactants, products, reversible, rateList, modificationRxn, originalRxnLine = ""):
        self.reactants = reactants[:]
        self.products = products[:]
        self.rateList = rateList[:]
        self.originalRxnLine = originalRxnLine

        # BNG, and consequently the rates in self.rateList use mass-action rate constants.
        # Therefore in the self.constructProperties function, we fill these in, using equation
        # (20) from "Modeling and Simulation of Intracellular Dynamics: ...".
        self.stochasticRateList = []
        

        ## Starting now, these will be generated based on the reactants,
        ## products, rateList information.

        self.reversible = False        
        self.decompositionRxn = False
        self.dimerizationRxn = False
        self.modificationRxn = False

        self.genericBinding = 0
        self.specificBinding = 0

        # This function sets correct values for the various properties.
        self.constructProperties()
        self.constructGenericBinding()

    def DEBUG_getReactionType(self):
        if self.decompositionRxn:
            return "DecompositionRxn"
        elif self.dimerizationRxn:
            return "DimerizationRxn"
        else:
            return "ModificationRxn"

    def resolveRates(self, aParameterInterpreter):
        for ndx in range(len(self.rateList)):
            if isinstance(self.rateList[ndx], type("hello world")):
                self.rateList[ndx] = aParameterInterpreter.evaluateExpression(self.rateList[ndx])

    def getGenericBinding(self):
        return self.genericBinding

    def constructGenericBinding(self):
        #### FORWARD
        # let rxn =
        # "Far1(MAPK_site,T306~none) + Fus3(target_site,T180~PO4,Y182~PO4) -> \Far1(MAPK_site!1,T306~none).Fus3(target_site!1,T180~PO4,Y182~PO4) kon_Fus3pTpY_Far1"
        # Then this function should return (("Far1", "MAPK_site"), ("Fus3, "target_site"))
        # This is the "generic binding that is formed.

        #### DECOMPOSITION
        # Let rxn =
        # "Far1(MAPK_site!1,T306~none).Fus3(target_site!1,T180~PO4,Y182~PO4) -> Far1(MAPK_site,T306~none) + Fus3(target_site,T180~PO4,Y182~PO4) koff_Fus3pTpY_Far1"
        # Then this function should return (("Far1", "MAPK_site"), ("Fus3, "target_site")).

        #### MODIFICATION
        # Let rxn = "Dig1(PO4_site~PO4) -> Dig1(PO4_site~none) 	kcat_nonspecific_dephosph"
        # Should return (Dig1, PO4_site, PO4), I *think*.
        

        if self.modificationRxn:
            self.genericBinding = util.calculateSignificantModificationSite(self.reactants, self.products)
            return self.genericBinding
        elif self.dimerizationRxn:
            self.genericBinding = util.calculateSignificantBinding(self.reactants, self.products)
        elif self.decompositionRxn:
            self.genericBinding = util.calculateSignificantBinding(self.products, self.reactants)
        else:
            raise ValueError

        self.genericBinding = list(self.genericBinding)
        self.genericBinding.sort()
        self.genericBinding = tuple(self.genericBinding)
        return self.genericBinding
    
    def constructSpecificBinding(self):
        #### FORWARD
        # let rxn =
        # "Far1(MAPK_site,T306~none) + Fus3(target_site,T180~PO4,Y182~PO4) -> \Far1(MAPK_site!1,T306~none).Fus3(target_site!1,T180~PO4,Y182~PO4) kon_Fus3pTpY_Far1"
        # This should return ("Far1(MAPK_site,T306~none)", 1, "MAPK_site"), ("Fus3(target_site,T180~PO4,Y182~PO4)", 1, "target_site")
        # This is the specific binding that is created.  

        #### DECOMPOSITION
        # Let rxn =
        # "Far1(MAPK_site!1,T306~none).Fus3(target_site!1,T180~PO4,Y182~PO4) -> Far1(MAPK_site,T306~none) + Fus3(target_site,T180~PO4,Y182~PO4) koff_Fus3pTpY_Far1"
        return
        
    def constructProperties(self):
        ## Set self.reversible
        if len(self.rateList) == 1:
            self.reversible = False
        elif len(self.rateList) == 2:
            self.reversible == True
        else:
            raise ValueError, "Error in constructing Reaction -- wrong length of rateList: " + str(len(self.rateList))

        # Because decomposition, dimerization, and modification are all
        # mutually exclusive, set them all in the same "if" block.

        if len(self.reactants) == 2 and len(self.products) == 1:
            self.dimerizationRxn = True
        elif len(self.reactants) == 1 and len(self.products) == 2:
            self.decompositionRxn = True
        elif len(self.reactants) == 1 and len(self.products) == 1:
            self.modificationRxn = True
        else:
            raise ValueError, "Error in constructing Reaction -- wrong length of reactants/products."

        self.calculateStochasticRateConstants()

    def setGlobalRxnVolume(RxnCls, newVolume): Reaction.volume = newVolume
    setGlobalRxnVolume = classmethod(setGlobalRxnVolume)

    def calculateStochasticRateConstants(self):
        # This is hacky and should be fixed.
        assert(Reaction.volume)

        ## Should this be self.volume?
        solutionVolume = Reaction.volume
        self.stochasticRateList.append( util.calculateStochasticRateConstantOnReactantArray(self.rateList[0], \
                                                                                            solutionVolume, \
                                                                                            self.getReactantStochiometryList()))

        if self.reversible:
            self.stochasticRateList.append( util.calculateStochasticRateConstantOnReactantArray(self.rateList[1], \
                                                                                                solutionVolume, \
                                                                                                self.getProductStochiometryList()))
    def getReactantStochiometryList(self):
        return self.getStochiometryFromList(self.reactants)

    def getProductStochiometryList(self):
        return self.getStochiometryFromList(self.products)

    def getStochiometryFromList(self, aList):
        theDict = {}

        for element in aList:
            theDict.setdefault(element, 0)
            theDict[element] += 1

        return theDict.values()
    
    def isModificationReaction(self):
        return self.modificationRxn

    def isDimerizationReaction(self):
        return self.dimerizationRxn

    def isDecompositionReaction(self):
        return self.decompositionRxn

    def debug(self):
        print self.originalRxnLine

        line0 = "Reactants: " + "  ".join(self.reactants)
        line1 = ""
        if self.reversible:
            line1 = "<->"
        else:
            line1 = "->"
        line2 = "Products: " + "  ".join(self.products)
        strRates = [str(x) for x in self.rateList]
        line3 = "Rates: " + " ".join(strRates)
        print line0
        print line1
        print line2
        print line3

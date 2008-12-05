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
# along with Moleculizer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Original Author:
#   Nathan Addy, Scientific Programmer	Voice: 510-981-8748
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#                     
#   
###############################################################################

import re, pdb, sys
from reaction import Reaction
from matchobject import MatchObject
from allosteryinterpreter import *
from namegenerator import NameGenerator
import util

class ReactionRules:
    def __init__(self, reactionRulesBlock, moleculeDictionary, parameterDictionary):
        self.rawReactionRules = reactionRulesBlock
        self.moleculeDictionary = moleculeDictionary
        self.parameterDictionary = parameterDictionary

        self.allosteryInterpreter = AllosteryInterpreter(self.moleculeDictionary)

        self.reactionRules = []
        self.initialize()


    def initialize(self):
        # First we set the static Volume for the Reaction class
        try:
            Reaction.setGlobalRxnVolume(self.parameterDictionary.getVolume())
        except:
            print "Error: the parameterDictionary cannot find a volume element"
            print "Crashing...."
            sys.exit(0)

        ## We classify all the reactions into a few classes for processing.
        isUnprocessable = lambda x: (x == 0) ## 0 is an unprocessable reaction
        nonProcessableReactions = [line for line in self.rawReactionRules if isUnprocessable(util.classifyReaction(line))]
        processibleReactions = [line for line in self.rawReactionRules if not isUnprocessable(util.classifyReaction(line))]

        ## NOTICE
        # There are currently 32 bad reactions, another 10 or so that were deleted at the outset.
        # ~ 8 Deletion reactions:  "Ste7 + Cell -> Cell kdeg_Ste7 DeleteMolecules"
        # ~ 1 Dilution Reaction:   "{MatchOnce}* + Cell -> Cell kdilution	exclude_reactants(1,Cell)"
        # ~ 23 Creation Reactions: "Cell -> Cell + Ste7(Ste5_site,MAPK_site,S359_T363~none) ksynth_Ste7"

        ## This takes each of the reactions in good reactions and expands them out.
        unidirectionalGoodReactions = []
        for line in processibleReactions:
            expandedLines = util.expandReactions(line)
            unidirectionalGoodReactions.extend(expandedLines)


        # THIS SHOULD ALL BE REDONE, because of DIMERIZATION + MODIFICATION RXNS

        # Get the dimerization reactions.
        tmpFunc = lambda x: (x == 1)
        dimerizationReactionsArray = [line for line in unidirectionalGoodReactions if tmpFunc(util.classifyReaction(line)) == True]

        # Get the decomposition reactions.
        tmpFunc = lambda x: (x == 2)
        decompositionReactionsArray = [line for line in unidirectionalGoodReactions if tmpFunc(util.classifyReaction(line)) == True]

        # Get the modification reactions.
        tmpFunc = lambda x: (x == 3)
        modificationReactionsArray = [line for line in unidirectionalGoodReactions if tmpFunc(util.classifyReaction(line)) == True]

        self.dimerizationReactions = {}
        self.decompositionReactions = {}
        self.modificationReactions = {}

        # Parse each reaction and install it in self.reactionRules.
        # Also determine whether it is a dimerization, decomposition, or modification reaction
        # and put it into the appropriate dictionary.
        for line in unidirectionalGoodReactions:
            aRxn = util.parseReaction(line, self.parameterDictionary)
            self.reactionRules.append(aRxn)
            
            gb = aRxn.getGenericBinding()
            if aRxn.isDimerizationReaction():
                if gb not in self.dimerizationReactions.keys():
                    self.dimerizationReactions[gb] = []
                self.dimerizationReactions[gb].append(aRxn)
            elif aRxn.isDecompositionReaction():
                if gb not in self.decompositionReactions.keys():
                    self.decompositionReactions[gb] = []
                self.decompositionReactions[gb].append(aRxn)
            elif aRxn.isModificationReaction():
                if gb not in self.modificationReactions.keys():
                    self.modificationReactions[gb] = []
                self.modificationReactions[gb].append(aRxn)

        for fullBindingCreated in self.dimerizationReactions:
            leftSite, rightSite = fullBindingCreated
            for rxn in self.dimerizationReactions[fullBindingCreated]:
                reactant1, reactant2 = rxn.reactants
                if reactant1.startswith(leftSite[0]):
                    self.allosteryInterpreter.createNewAllosteryAtSite(leftSite, reactant1)
                    self.allosteryInterpreter.createNewAllosteryAtSite(rightSite, reactant2)
                elif reactant1.startswith(rightSite[0]):
                    self.allosteryInterpreter.createNewAllosteryAtSite(rightSite, reactant1)                    
                    self.allosteryInterpreter.createNewAllosteryAtSite(leftSite, reactant2)
                else:
                    raise ValueError

        for fullBindingCreated in self.decompositionReactions.keys():
            leftSite, rightSite = fullBindingCreated
            for rxn in self.decompositionReactions[fullBindingCreated]:
                product1, product2 = rxn.products
                if product1.startswith(leftSite[0]):
                    self.allosteryInterpreter.createNewAllosteryAtSite(leftSite, product1)
                    self.allosteryInterpreter.createNewAllosteryAtSite(rightSite, product2)
                elif product1.startswith(rightSite[0]):
                    self.allosteryInterpreter.createNewAllosteryAtSite(rightSite, product1)                    
                    self.allosteryInterpreter.createNewAllosteryAtSite(leftSite, product2)
                else:
                    raise ValueError

        self.writeAllostericSitesToMoleculeObject()

    def writeAllostericSitesToMoleculeObject(self):
        for bindingSite in self.allosteryInterpreter.allosteryMapping.keys():
            containingMol, bindingSiteName = bindingSite
            for inducingComplex in self.allosteryInterpreter.allosteryMapping[bindingSite]:
                siteShapeName = self.allosteryInterpreter.allosteryMapping[bindingSite][inducingComplex]
                # Now lookup that mol in self.moleculeDictionary.registeredMoleculesDictionary
                try:
                    if not siteShapeName in self.moleculeDictionary.registeredMoleculesDictionary[containingMol].bindingSites[bindingSiteName]:
                        self.moleculeDictionary.registeredMoleculesDictionary[containingMol].bindingSites[bindingSiteName].append(siteShapeName)
                except:
                    print "ERRORERRORERROR"
                    pdb.set_trace()
                    a = 1
                      
    def writeReactionGensElement(self, reactionGensElement):
        allAvailableBindingReactions = []
        length = len(self.dimerizationReactions.keys())
        allAvailableBindingReactions.extend(self.dimerizationReactions.keys())
        allAvailableBindingReactions.extend(self.decompositionReactions.keys())
        allAvailableBindingReactions = list(set(allAvailableBindingReactions))

        for bindingType in allAvailableBindingReactions:
            dimerizationGenElmt = XmlObject("dimerization-gen")
            dimerizationGenElmt.attachToParent(reactionGensElement)
            self.writeDimerizationGen(bindingType, dimerizationGenElmt)

        for modificationSite in self.modificationReactions:
            for rxn in self.modificationReactions[modificationSite]:
                self.addModificationRxnToReactionGenElement(rxn, reactionGensElement, modificationSite)
        return


    def addModificationRxnToReactionGenElement(self, rxn, reactionGensElement, modificationInfo):

        # Create all the objects
        omniGenElmt = XmlObject("omni-gen")
        enablingOmniPlexElmt = XmlObject("enabling-omniplex")
        smallMolExchangesElmt = XmlObject("small-mol-exchanges")
        modificationExchangeElmt = XmlObject("modification-exchanges")
        additionalReactantSpeciesElmt = XmlObject("additional-reactant-species")
        additionalProductSpeciesElmt = XmlObject("additional-product-species")
        rateElmt = XmlObject("rate")

        # Now add them all to the xml structure

        omniGenElmt.attachToParent(reactionGensElement)
        enablingOmniPlexElmt.attachToParent(omniGenElmt)
        smallMolExchangesElmt.attachToParent(omniGenElmt)
        modificationExchangeElmt.attachToParent(omniGenElmt)
        #additionalReactantSpeciesElmt.attachToParent(omniGenElmt)
        #additionalProductSpeciesElmt.attachToParent(omniGenElmt)

        # Add a plex element to the enablingOmniPlexElmt
        plexElmt = XmlObject("plex")
        plexElmt.attachToParent(enablingOmniPlexElmt)
        assert(len(rxn.reactants) == 1 and len(rxn.products) == 1)
        reactant = rxn.reactants[0]
        reactantArray = reactant.split(".")

        molNdxToInstanceName = {}
        molNdxToRefName = {}
        tmpNameGenerator = {}
        bindingSites = {}
        unboundSites = []
        modificationSites = []
        for ndx in range(len(reactantArray)):
            molSpec = reactantArray[ndx]
            parenNdx = molSpec.find("(")
            molRefName = molSpec[:parenNdx]
            if not molRefName in tmpNameGenerator.keys():
                molInstanceName = "the-" + molRefName
                tmpNameGenerator[molRefName] = 1
            else:
                molInstanceName = molRefName + "_" + str(tmpNameGenerator[molRefName])
                tmpNameGenerator[molRefName] += 1

            molNdxToRefName[ndx] = molRefName
            molNdxToInstanceName[ndx] = molInstanceName

            # Create the mol-intance element
            molInstanceElmt = XmlObject("mol-instance")
            molInstanceElmt.addAttribute("name", molInstanceName)
            molRefElmt = XmlObject("mol-ref")
            molRefElmt.addAttribute("name", molRefName)
            molRefElmt.attachToParent(molInstanceElmt).attachToParent(plexElmt)
            
            # Now delve into the bindings of this muthafucka.
            # Note bound components, unbound componends, and modifications...

            bindingInformation = molSpec[parenNdx + 1: -1]
            bindingInfoArray = bindingInformation.split(",")

            bindings = [x for x in bindingInfoArray if "~" not in x]
            modifications = [x for x in bindingInfoArray if "~" in x]
            boundBindings = [x for x in bindings if "!" in x and x[-1] != "+"]
            unboundBindings = [x for x in bindings if "!" not in x]

            # Process modifications
            for modification in modifications:
                site, state = modification.split("~")
                modificationSites.append( (ndx, site, state))

            # Process unbound binding sites
            for unboundBinding in unboundBindings:
                unboundSites.append( (ndx, unboundBinding))

            for boundHalfBinding in boundBindings:
                site, bindingNdx = boundHalfBinding.split("!")
                if not bindingNdx in bindingSites.keys():
                    bindingSites[bindingNdx] = []
                bindingSites[bindingNdx].append( (ndx, site))
                
        # End initial pass, now we can put in all the bindingSites,
        # unbound bindingSites, and modification states.

        for ndx in bindingSites.keys():
            bindingElmt = XmlObject("binding")
            bindingElmt.attachToParent(plexElmt)
            for halfBinding in bindingSites[ndx]:
                molInstanceRefElmt = XmlObject("mol-instance-ref")
                molInstanceRefElmt.addAttribute("name", molNdxToInstanceName[halfBinding[0]])

                bindingSiteRefElmt = XmlObject("binding-site-ref")
                bindingSiteRefElmt.addAttribute("name", halfBinding[1])

                bindingSiteRefElmt.attachToParent(molInstanceRefElmt).attachToParent(bindingElmt)

        # I think I should try to filter out the modification(the relevant one)
        if len(modificationSites) > 0:
            instanceStatesElmt = XmlObject("instance-states")
            instanceStatesElmt.attachToParent(plexElmt)
            
        for modMolNdx in set([x[0] for x in modificationSites]):
            modMap = [ (x[1], x[2]) for x in modificationSites if x[0] == modMolNdx]

            modMolInstanceRefElmt = XmlObject("mod-mol-instance-ref")
            modMolInstanceRefElmt.addAttribute("name", molNdxToInstanceName[modMolNdx])

            modMapElmt = XmlObject("mod-map")
            modMapElmt.attachToParent(modMolInstanceRefElmt).attachToParent(instanceStatesElmt)

            for modification in modMap:
                modSiteRefElmt = XmlObject("mod-site-ref")
                modSiteRefElmt.addAttribute("name", modification[0])

                modRefElmt = XmlObject("mod-ref")
                modRefElmt.addAttribute("name", modification[1])
                modRefElmt.attachToParent(modSiteRefElmt).attachToParent(modMapElmt)

        if len(unboundSites) > 0:
            unboundSitesElmt = XmlObject("unbound-sites")
            unboundSitesElmt.attachToParent(plexElmt)
        for unboundSite in unboundSites:
            instanceRefElmt = XmlObject("instance-ref")
            instanceRefElmt.addAttribute("name", molNdxToInstanceName[unboundSite[0]])

            siteRefElmt = XmlObject("site-ref")
            siteRefElmt.addAttribute("name", unboundSite[1])
            siteRefElmt.attachToParent(instanceRefElmt).attachToParent(unboundSitesElmt)


        ## Now create the modification exchange
        molRef, modSite, modInitialState = modificationInfo
        realNdx = -1
        for ndx in range(len(reactantArray)):
            if reactantArray[ndx].startswith(molRef) and (modSite + "~" + modInitialState) in reactantArray[ndx]:
                realNdx = ndx
                break
        else:
            pdb.set_trace()
            a = 1
            raise ValueError
        
        modMolInstanceRefElmt = XmlObject("mod-mol-instance-ref")
        modMolInstanceRefElmt.addAttribute("name", molNdxToInstanceName[realNdx])

        modSiteRefElmt = XmlObject("mod-site-ref")
        modSiteRefElmt.addAttribute("name", modSite)

        modSiteRefElmt.attachToParent(modMolInstanceRefElmt).attachToParent(modificationExchangeElmt)

        installedModRefElmt = XmlObject("installed-mod-ref")
        # To find the result, we must find out what state the modSite modification is in the realNdx-th plex on the product side
        modifications = util.getModificationList(rxn.products[0].split('.')[realNdx])
        endState = ""
        for x in modifications:
            if x.startswith(modSite):
                endState = x.split('~')[1]
                break
        else:
            raise ValueError

        installedModRefElmt.addAttribute("name", endState)
        installedModRefElmt.attachToParent(modificationExchangeElmt)
        rateElmt.addAttribute("value", self.parameterDictionary.evaluateExpression(rxn.rateList[0]))
        rateElmt.attachToParent(omniGenElmt)

    def writeDimerizationGen(self, bindingType, dimerizationGenElmt):
        for halfBinding in bindingType:
            molRefElmt = XmlObject("mol-ref")
            molRefElmt.addAttribute("name", halfBinding[0])
            siteRefElmt = XmlObject("site-ref")
            siteRefElmt.addAttribute("name", halfBinding[1])
            siteRefElmt.attachToParent(molRefElmt).attachToParent(dimerizationGenElmt)
            
        defaultOnRateElmt = XmlObject("default-on-rate")
        defaultOffRateElmt = XmlObject("default-off-rate")
        defaultOnRateElmt.attachToParent(dimerizationGenElmt)
        defaultOffRateElmt.attachToParent(dimerizationGenElmt)

        defaultCreated = False

        # Maybe here I should write a function to process all the shit
        # for this bindingType.
        # It could return a series of 4-tuples:
        # (binding_allostery_on_left, binding_allostery_on_right,\
        #   on_rate, off_rate)
        rxnsOnBinding = self.constructReactionAllosteryRates(bindingType)
        invVolume = 1.0 / self.parameterDictionary.getVolume()
        for allostericReactionOnBinding in rxnsOnBinding:
            left_allostery, right_allostery, onRate, offRate = allostericReactionOnBinding
            if left_allostery == "default" and right_allostery == "default":
                defaultOnRateElmt.addAttribute("value", invVolume * self.parameterDictionary.evaluateExpression(onRate))
                defaultOffRateElmt.addAttribute("value", invVolume * self.parameterDictionary.evaluateExpression(offRate))
                defaultCreated = True
            else:
                # Process dimerization rxn
                alloRatesElmt = XmlObject("allo-rates")
                alloRatesElmt.attachToParent(dimerizationGenElmt)
                
                siteShapeRef1Elmt = XmlObject("site-shape-ref")
                siteShapeRef1Elmt.addAttribute("name", left_allostery)
                siteShapeRef1Elmt.attachToParent(alloRatesElmt)

                siteShapeRef2Elmt = XmlObject("site-shape-ref")
                siteShapeRef2Elmt.addAttribute("name", right_allostery)
                siteShapeRef2Elmt.attachToParent(alloRatesElmt)

                onRateElmt = XmlObject("on-rate")
                onRateElmt.addAttribute("value", invVolume * self.parameterDictionary.evaluateExpression(onRate))
                onRateElmt.attachToParent(alloRatesElmt)
                
                offRateElmt = XmlObject("off-rate")
                offRateElmt.addAttribute("value", invVolume * self.parameterDictionary.evaluateExpression(offRate))
                offRateElmt.attachToParent(alloRatesElmt)
        else:
            if not defaultCreated:
                defaultOnRateElmt.addAttribute("value", 0.0)
                defaultOffRateElmt.addAttribute("value", 0.0)
            

##             if not defaultCreated:
##                 pdb.set_trace()
##                 a = 1
##                 raise ValueError, "Error, no reaction rates were set for this class of dimerization reactions."



    def getSaturationRateConstant(self, saturationExpressionToken):
        ## This needs to be implemented
        return 1.0

    def addAllostericOmnis(self, allostericOmnisElement):
        self.allosteryInterpreter.addAllostericOmnis(allostericOmnisElement)

    def constructReactionAllosteryRates(self, bindingType):
        # This function should look up bindingType in both self.dimerizationReactions
        # (reactions that construct that binding) and self.decompositionReactions
        # (reactions that break up that binding), and should try to pair them.

        leftHalfBinding, rightHalfBinding = bindingType

        rxnsDimerizingBindingType = self.dimerizationReactions[bindingType]
        rxnsDecomposingBindingType = self.decompositionReactions[bindingType]
        


##         if len(rxnsDimerizingBindingType) == 4:
##             pdb.set_trace()


        pairedRxns = []        
        usedNdxs = []
        for dimRxn in rxnsDimerizingBindingType:
            reactants = dimRxn.reactants
            for ndx in range(len(rxnsDecomposingBindingType)):
                ## This is bad methinks.  Probably not, but perhaps.
                if ndx in usedNdxs:
                    continue
                products = rxnsDecomposingBindingType[ndx].products
                if (reactants[0] == products[0] and reactants[1] == products[1]) or (reactants[0] == products[1] and reactants[1] == products[0]):
                    usedNdxs.append(ndx)
                    pairedRxns.append( (dimRxn, rxnsDecomposingBindingType[ndx]))
                    break
            else:
                # No match for the dimerization reaction was found.  For now just print out a "barf"
                raise ValueError, "Unmatched reaction..."


        finalStuff = []
        for reactionPair in pairedRxns:
            assert(not reactionPair[0].reversible and not reactionPair[1].reversible)
            
            allostericComplexes = reactionPair[0].reactants
            alloPlex1, alloPlex2 = allostericComplexes
            if util.matches(leftHalfBinding, alloPlex1) and util.matches(rightHalfBinding, alloPlex2):
                leftAllostericShape = self.allosteryInterpreter.allosteryMapping[leftHalfBinding][alloPlex1]
                rightAllostericShape = self.allosteryInterpreter.allosteryMapping[rightHalfBinding][alloPlex2]
                onRate = reactionPair[0].stochasticRateList[0]
                offRate = reactionPair[1].stochasticRateList[0]
                point = (leftAllostericShape, rightAllostericShape, onRate, offRate)
                finalStuff.append(point)
            elif util.matches(leftHalfBinding, alloPlex2) and util.matches(rightHalfBinding, alloPlex1):
                leftAllostericShape = self.allosteryInterpreter.allosteryMapping[leftHalfBinding][alloPlex2]
                rightAllostericShape = self.allosteryInterpreter.allosteryMapping[rightHalfBinding][alloPlex1]
                onRate = reactionPair[0].stochasticRateList[0]
                offRate = reactionPair[1].stochasticRateList[0]
                point = (leftAllostericShape, rightAllostericShape, onRate, offRate)
                finalStuff.append(point)
            else:
                pdb.set_trace()
                a = 1
                raise ValueError
        return finalStuff

            
            
        



                                                                                                     
            
        

        

        
        
        

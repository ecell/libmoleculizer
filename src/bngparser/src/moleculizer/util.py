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

from reaction import *
import pdb

N_Avagadro = 6.0221415e23

def sort(sequence, func = 0):
    if func:
        sequence.sort(func)
    else:
        sequence.sort()
    return sequence

def prod(aList):
    product = 1
    for num in aList:
        product *= num


def matches(halfBinding, complex):
    if halfBinding[0] in complex and halfBinding[1] in complex:
        return True
    else:
        return False

def calculateStochasticRateConstantOnReactantArray(gmaRateConstant, solutionVolume, stochiometryList):
    # This conversion was found in Wolkenhauer et al --: Modeling and Simulation of Intracellular Dynamics:
    # Choosing an Appropriate Framework", equation 20.

    molecularityMinusOne = sum(stochiometryList) - 1
    denom = (N_Avagadro * solutionVolume) ** molecularityMinusOne
    leftHandTerm = gmaRateConstant / denom

    mult = lambda x, y: x * y
    fact = lambda x: reduce(mult, range(1, x + 1))
    rightHandTerm = reduce(mult, map(fact, stochiometryList))

    return leftHandTerm * rightHandTerm 

def expandReactions(rxnRuleLine):
    if not "<->" in rxnRuleLine:
        return [rxnRuleLine]

    rxn = rxnRuleLine.split()
    
    assert("Cell" not in rxn)
    assert("{MatchOnce}*" not in rxn)
    assert("DeleteMolecules" not in rxn)
    
    rates = []
    rates.append(rxn.pop())
    rates.append(rxn.pop())
    rates.reverse()

    midNdx = rxn.index("<->")

    reactants = rxn[:midNdx]
    reactants = [x for x in reactants if x != "+"]

    products = rxn[midNdx + 1:]
    products = [x for x in products if x != "+"]

    forwardReactionLine = " -> ".join([ " + ".join(reactants), " + ".join(products)]) + " " + str(rates[0])
    reverseReactionLine = " -> ".join([ " + ".join(products), " + ".join(reactants)]) + " " + str(rates[1])

    return [forwardReactionLine, reverseReactionLine]
    

def parseReaction(rxnRuleLine, paramInterpreter = None):
    ## This should be totally refactored out into the Reaction function.  Yuck.
    
    ## This is given a reaction with classification (see function "classifyReaction"
    ## that is ge 1 and le 6.)

    reversible = -1
    reactants = -1
    products = -1
    rateList = []
    modificationRxn = -1
    ndx = -1

    rxnComponents = rxnRuleLine.strip().split()
    if "<->" in rxnComponents:
        reversible = True
        rateList.append(rxnComponents.pop())
        rateList.append(rxnComponents.pop())
        rateList.reverse()
        ndx = rxnComponents.index("<->")
    elif "->" in rxnComponents:
        reversible = True
        rateList.append(rxnComponents.pop())
        ndx = rxnComponents.index("->")
    else:
        raise ValueError, ("Error, reaction line apparently has no reaction symbols...\n\t" + rxnRuleLine)

    reactants = rxnComponents[:ndx]
    reactants = [x for x in reactants if x != "+"]
    products = rxnComponents[ndx + 1:]
    products = [x for x in products if x != "+"]

    if len(reactants) == len(products) == 1:
        modificationRxn = True
    else:
        modificationRxn = False

    if paramInterpreter:
        for ndx in range(len(rateList)):
            if isinstance(rateList[ndx], type("string type")):
                rateList[ndx] = paramInterpreter.evaluateExpression(rateList[ndx])

    ## TODOTesting: Ensure that the values going into this constructor are all correct.
    # pdb.set_trace()
    return Reaction(reactants, products, reversible, rateList, modificationRxn, rxnRuleLine)

def classifyReaction(rxnRuleLine):
    ### 0: bad reaction
    ### 1: 2->1 rxn
    ### 2: 1->2 rxn
    ### 3: modificationrxn


    ### IMPORTANT!!!
    ## It may seem tempting to rewrite this function in terms of
    ## the parseReaction function above.  BEWARE -- Dragons lurk
    ## there.  The parseReaction function assumes a reaction line
    ## such that 1 <= classifyReaction(line) <= 6.
    

    if "Cell" in rxnRuleLine or "DeleteMolecules" in rxnRuleLine or "{MatchOnce}*" in rxnRuleLine:
        return 0
    
    rxn = rxnRuleLine.split()
    midNdx = -1
    if "<->" in rxn:
        rxn = rxn[:-2]
        midNdx = rxn.index("<->")
    elif "->" in rxn:
        rxn = rxn[:-1]
        midNdx = rxn.index("->")
    else:
        raise ValueError, "No rxn symbol is found in this rxn rule line...\n\t" + rxnRuleLine
    
    reactants = rxn[:midNdx]
    reactants = [x for x in reactants if x != "+"]
    
    products = rxn[midNdx + 1:]
    products = [x for x in products if x != "+"]

    if len(reactants) == 1 and len(products) == 2:
        return 2
    elif len(reactants) == 2 and len(products) == 1:
        return 1
    elif len(reactants) == 1 and len(products) == 1:
        return 3

def getBindingSites(molSpec):
    molSpec.strip()

    bndSiteNdx = molSpec.find('(')
    if bndSiteNdx == -1:
        return []
    
    molSpec = molSpec[bndSiteNdx + 1:]
    molSpec = molSpec[:-1]

    bngMolComponents = molSpec.split(',')
    bngMolComponents = [x for x in bngMolComponents if "~" not in x]
    return bngMolComponents

def bindingSpecificationRepresentsBound(bindingComponent):
    return "!" in bindingComponent

def getBindingNdxFromBoundComponent(bindingComponent):
    ndx = bindingComponent.index('!')
    return bindingComponent[ndx + 1:].strip()

def getMolNameFromFullMolSpec( bngMolSpecification ):
    # This function takes a MOL name - Ste11(bnd_1, bnd_2, a~0~1) and returns
    # the name, in this case Ste11
    
    bngMolSpecification = bngMolSpecification.strip()

    ndx = bngMolSpecification.find('(')
    if ndx == -1:
        return bngMolSpecification
    else:
        return bngMolSpecification[:ndx]


def getBindingSpecificationFromNdx(bngComplexSpecification, bindingNdx):
    # This function takes a complex name like Ste5(ste_11!1, ste_12, ste13).Ste11(ste_5!1, b~0~1)
    # and an int like 1 and reurns [(Ste5, ste11), (Ste11, ste_5)].
    
    bindingNdxStr = str(bindingNdxStr)
    
    bngMolComponentList = bngComplexSpecification.split('.')
    bindingComponents = []
    
    for component in bngMolecularComponentList:
        if bindingNdxStr in component:
            molName = getMolNameFromMolSpec( component )
            bindingName = getBindingNameWithBoundID( component, bindingNdxStr )
            bindingComponents.append( (molName, bindingName) )
    else:
        assert( len(bindingComponents) == 2)
        bindingComponents.sort()
        return bindingComponents
    raise ValueError, "Error, bad binding.  Plus I need a better error message..."


def calculateSignificantBinding(reactantList1, reactantList2):

    smallerList, longerList = sort( [reactantList1, reactantList2], lambda a, b : cmp(len(a), len(b)))

    ## Collate up the binding numbers from the side with fewer bindings.
    bindingNumList = []
    for participant in longerList:
        for bindingSite in getBindingSites(participant):
            if bindingSpecificationRepresentsBound(bindingSite):
                bindingNdx = getBindingNdxFromBoundComponent(bindingSite)

                if not bindingNdx in bindingNumList:
                    bindingNumList.append(bindingNdx)

    missingNdx = 0
    bindingCreated = 0

    # Now cycle through the binding sites of the larger complex, to find the newly
    # created binding.
    for mol in smallerList[0].split('.'):
        for bindingSite in getBindingSites(mol):
            if bindingSpecificationRepresentsBound(bindingSite) and getBindingNdxFromBoundComponent(bindingSite) not in bindingNumList:
                missingNdx = getBindingNdxFromBoundComponent(bindingSite)
                bindingCreated = getSpecificBindingRepresentation(smallerList[0], missingNdx)

                # This may now have length one or 2
                bindingCreated = list(bindingCreated)
                bindingCreated.sort()
                return tuple(bindingCreated)


    # This binding type probably just is one for which no binding is created....
    if "Cell" in reactantList1:
        return []
    else:
        pdb.set_trace()
        a = 1
        raise ValueError


def calculateSignificantModificationSite(reactantList1, reactantList2):
    assert( len(reactantList1) == len(reactantList2) == 1 )

    pre = reactantList1[0]
    preList = pre.split('.')

    post = reactantList2[0]
    postList = post.split('.')

    modificationDict = {}
    modifiedModificationsArray = []

    for ndx in range(len(preList)):
        modificationDict[ndx] = {}
        for value in util.getModificationList(preList[ndx]):
            modificationSite, modificationState = value.split("~")
            modificationDict[ndx][modificationSite] = modificationState

    for ndx in range(len(postList)):
        for value  in util.getModificationList(postList[ndx]):
            modificationSite, modificationState = value.split("~")
            if modificationDict[ndx][modificationSite] != modificationState:
                modifiedModificationsArray.append( (util.getMolNameFromFullMolSpec(postList[ndx]), modificationSite, modificationDict[ndx][modificationSite]))

    assert(len(modifiedModificationsArray) == 1)
    return modifiedModificationsArray[0]



def representsNullBinding(bindingName):
    # Is this function used? 
    listOfNullNames = ["none"]
    return bindingName.lower() in listOfNullNames

def getSpecificBindingRepresentation(complexSpecification, index):
    bindingRepresentation = []
    
    componentList = complexSpecification.split('.')
    for component in componentList:
        if ("!" + index ) in component:
            molName = getMolNameFromFullMolSpec(component)
            bindingName = getBindingNameByIndex(component, index)
            bindingRepresentation.append( (molName, bindingName) )

    bindingRepresentation.sort()
    return tuple(bindingRepresentation)

def getBindingNameByIndex(component, index):
    ndx = component.index('(')
    derivedComponent = component[ndx + 1:]
    derivedComponent = derivedComponent[:-1]
    bindingsList = derivedComponent.split(',')
    for binding in bindingsList:
        if ("!" + index) in binding:
            try:
                excNdx = binding.index('!')
                return binding[:excNdx]
            except:
                pdb.set_trace()
                x = getBindingNameByIndex(component, index)
                return x
    else:
        pdb.set_trace()
        a = 1
        raise ValueError
    
def getModification(modificationRepresentation):
    modificationRepresentation = modificationRepresentation.strip()
    assert('~' in modificationRepresentation)

    modificationSite, modificationState = modificationRepresentation.split('~')
    return modificationSite, modificationState

def getModificationList(molSpecification):
    ndx = molSpecification.find('(')
    if ndx == -1:
        return []

    molSpecification = molSpecification[ndx + 1:]
    molSpecification = molSpecification[:-1]
    bindmodList = molSpecification.split(',')
    modificationList = [x for x in bindmodList if "~" in x]
    return modificationList
        

def getUniqueList(aList):
    return list(set(aList))

def getModificationBindingSites(reactant, product):
    return []

def getAllostericBindingSites(longList, smallList):
    # Need to return two pairs, ((bindingSiteName, Mol), Complex), ((bindingSiteName, Mol), Complex)

    if not (len(longList) == 2 and len(smallList) == 1):
        raise ValueError, "Error, reactants and products don't have the right length..."
    
    allosterySpecifications = []

    aBinding = calculateSignificantBinding(smallList, longList)
    for halfBinding in aBinding:
        molName, bindingSite = halfBinding
        # name should exist in exactly one of the reactant complexes
        # this is unrealistic but ok for a first approximation (it works within
        # the constraints of *THE* model we've been given.

        array = [plex.startswith(molName) for plex in longList]

## Comment this out for the time being
##         try:
##             assert(array.count(True) == 1)
##         except:
##             pdb.set_trace()
##             a = 10

        theComponent = 0
        for component in longList:
            if molName in component:
                theComponent = component
                break
        else:
            raise ValueError, "Component not found..."

        allosterySpecifications.append( ((bindingSite, molName), theComponent))

    return allosterySpecifications


def isExcludedReaction(aLine):
    components = aLine.split()
    if "Cell" in components:
        return True
    else:
        return False

def compare(plex, moleculeDictionary):
    if '.' in plex:
        return False

    if "!" in plex:
        return False

    if plex.count("~") != plex.count("~none"):
        return False

    return True

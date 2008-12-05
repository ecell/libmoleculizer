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

import random
import pdb

class BngMoleculeDefinition:
    class BadMolDefinitionException(Exception): pass
    def __init__(self, bngMoleculeDefinitionLine):
	self.name = ""
	self.bindingSites = []
	self.modificationSites = []
	self.modificationStates = []
	self.parseBNGMoleculeDefinition(bngMoleculeDefinitionLine)
 
    def parseBNGMoleculeDefinition(self, moleculeDefinition):
	parenIndex = moleculeDefinition.find('(')
	if parenIndex == -1: # The definition is only a name.
	    self.name = moleculeDefinition.strip()
	    return
	self.name, unparsedBindings = moleculeDefinition[:parenIndex], moleculeDefinition[parenIndex:]
	unparsedBindings = unparsedBindings[1:]
	unparsedBindings = unparsedBindings[:-1]
	unparsedBindingsArray = unparsedBindings.split(',')
	bindingSites = [x for x in unparsedBindings if '~' not in x]
	modificationSites = [x for x in unparsedBindings if '~' in x]
	for bindingSite in bindingSites:
	    self.addBindingSite(bindingSite)
	for modificationSite in modificationSites:
	    self.addModificationSite(modificationSite)
	    
    def addModificationSite(self, modificationSite):
	# Modification site will be in form: name~state1~state2[~staten]*n
	parsedModificationSite = modificationSite.split('~')
	modificationSiteName = parsedModificationSite.pop(0)
	for name, default in self.modificationSites:
	    if modificationSiteName == name:
		raise BngMoleculeDefinition.BadMolDefinitionException("Error: Modification Site '%s' defined twice in molecule '%s'" % (modificationSite, self.getName()))
	self.modificationSites.append(modificationSiteName, parsedModificationSite[0])
	modificationStates = [state for state in parsedModificationSite if state not in self.modificationStates]
	self.modificationStates.extend(modificationStates)
    
    def addBindingSite(self, bindingSite):
	if bindingSite in self.bindingSites:
	    raise BngMoleculeDefinition.BadMolDefinitionException("Error: Binding Site '%s' defined twice in molecule '%s'" % (bindingSite, self.getName()))
	self.bindingSites.append(bindingSite)
	
    def getName(self):
	return self.name

    def getBindingSites(self):
	return self.bindingSites[:]

    def getModificationSites(self):
	return self.modificationSites[:]

    def getModificationStates(self):
	return self.modificationStates[:]

    def isLegalBinding(self, bindingName):
	return (bindingName in self.bindingSites)
    def isLegalModification(self, modification):
	# Modification can be specified as (modsite, value) or as modsite~value
	if isinstance(modification, str):
	    modSite, modState = modification.split('~')
	elif isinstance(modification, tuple) and len(modification)==2:
	    modSite, modState = modification
	else:
	    raise BngMoleculeDefinition.BadMolDefinitionException("Error: Modification '%s' is not a valid specification of a modification" % modification)
	isLegal = False
	for definedModification, default in self.getModificationSites():
	    if definedModification == modification:
		isLegal == True
	if isLegal and modState in self.modificationStates:
	    return True
	else:
	    return False

class BngSpeciesInstantiation:
    class Namer:
        def __init__(self):
            self.names = {}
        def getUniqueInstance(self, molReference):
            if not self.names.has_key(molReference): self.names[molReference] = 0
            if self.names[molReference] == 0:
                uniqueInstanceName = "the-" + molReference
            else:
                uniqueInstanceName = molReference + "-" + str(self.names[molReference])
            self.names[molReference] += 1
            return uniqueInstanceName
            
    def __init__(self, bngSpeciesLine):
        self.bngDefinition = bngSpeciesLine

	self.moleculeRefsInComplex = {}
        self.moleculeInstancesInComplex = {}

	self.bindings = dict()
	self.modificationsByMoleculeIndex = []
        
	self.regexBinding = []

	self.parseBNGSpeciesDefinition(bngSpeciesLine)


    def getBNGSpecification(self):
        return self.bngDefinition

    def getSize(self):
        return len(self.moleculeRefsInComplex.keys())

    def constructEntireComplexInstanceName(self):
        if self.getSize() == 1:
            return self.moleculeRefsInComplex[self.moleculeRefsInComplex.keys()[0]] + "-singleton"

        instanceName = ""
        for ndx in self.moleculeRefsInComplex:
            instanceName += (self.moleculeRefsInComplex[ndx] + "_")
        else:
            instanceName = instanceName[:-1] + "-complex-" + str(random.randint(1, 1e6))
            return instanceName
	
    def parseBNGSpeciesDefinition(self, bngSpeciesLine):
        __namerObj = BngSpeciesInstantiation.Namer()

	complexedMolecules = bngSpeciesLine.split('.')

	for moleculeNdx in range(len(complexedMolecules)):
	    molecule = complexedMolecules[moleculeNdx]

	    # parse the molecule for bindings and modifications
	    if '(' not in molecule:
		self.moleculeRefsInComplex[moleculeNdx] = molecule
                self.moleculeInstancesInComplex[moleculeNdx] = __namerObj.getUniqueInstance(molecule)
		continue

 	    parenNdx = molecule.find('(')
            moleculeTypeReference = molecule[:parenNdx]
            moleculeTypeInstance = __namerObj.getUniqueInstance(moleculeTypeReference)
            stateInfo = molecule[parenNdx+1:]
	    stateInfo = stateInfo[:-1]
	    stateInfo = stateInfo.split(',')
            bindings = [x for x in stateInfo if "~" not in x]
	    modifications = [x for x in stateInfo if "~" in x]

	    self.moleculeRefsInComplex[moleculeNdx] = moleculeTypeReference
            self.moleculeInstancesInComplex[moleculeNdx] = moleculeTypeInstance

	    for mod in modifications:
		# Each should be in the form "mod~val"
		mod, value = mod.split("~")
		self.modificationsByMoleculeIndex.append((moleculeNdx, mod, value))

	    for binding in bindings:
		if "!" in binding:
		    bindingInfo = binding.split("!")
		    assert( len(bindingInfo) == 2)

		    isint = False
		    try:
			int(bindingInfo[1])
			isint = True
		    except:
			isint = False
	
		    if isint:
			bindingNumber = int(bindingInfo[1])
			bindingInfo = (moleculeNdx, bindingInfo[0])
			if bindingNumber in self.bindings.keys():
			    self.bindings[bindingNumber].append( bindingInfo )
			else:
			    self.bindings[bindingNumber] = [ bindingInfo ]
		    elif bindingInfo[1] == '+':
			# Handle this strange case
			# maybe add it to a MandatoryUnspecifiedBinding?
			bindingInfo = (moleculeNdx, bindingInfo[0])
			self.regexBinding.append(bindingInfo)
		    else:
			print "Unhandled thing"
			print bindingInfo[1]
			raise "unhandled binding type"
			

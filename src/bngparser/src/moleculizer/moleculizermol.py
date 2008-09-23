###############################################################################
# BNGMZRConverter - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007, 2008 The Molecular Sciences Institute
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
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
class MoleculizerMol:

    class BadMolDefinitionException(Exception): pass
    
    modificationStates = []
    
    def __init__(self, aLine):
	self.name = ""
	self.bindingSites = {}
	self.modificationSites = []
	self.parseMolLine(aLine)

    def setName(self, aName):
	self.name = aName

    def getName(self):
	if self.name.strip() == "":
	    raise MoleculizerMol.BadMolDefinitionException("Error: Mol has no name")
	else:
	    return self.name
    def addBindingSite(self, bindingSite):
        ## A binding site consists only of its name -- ie "Ste11_site" or some such.
	if not bindingSite in self.bindingSites:
            # Keep an eye on this line, not positive what I intended this to mean....
	    self.bindingSites[bindingSite] = ["default"]

	else:
	    raise MoleculizerMol.BadMolDefinitionException("Error: Mol '%s' already has the binding site '%s'" % (bindingSite))

    def addBindingSiteShape(self, bindingSite, bindingShape):
	if bindingSite not in self.bindingSites.keys():
	    raise MoleculizerMol.BadMolDefinitionException("Error: Mol '%s' has no binding site named '%s'" % (self.name, bindingSite))
	if bindingShape in self.bindingSites[bindingSite]:
	    raise MoleculizerMol.BadMolDefinitionException("Error: Mol '%s' already has a binding shape named '%s' at binding site '%s'" % (self.name, bindingSite, bindingShape))
	self.bindingSites[bindingSite].append(bindingShape)

    def addModificationSite(self, modification):
	# Takes in a string like "T034~U~P"; the first is taken as the default.

	modification = modification.strip()
	modificationArray = modification.split('~')
	modificationName = modificationArray.pop(0)
	if (modificationName, modificationArray[0]) in self.modificationSites:
	    raise MoleculizerMol.BadMolDefinitionException("Error: Mol '%s' already has a modification site '%s'" % (self.name, modificationName))
	for mod in modificationArray:
	    if mod not in MoleculizerMol.modificationStates:
		MoleculizerMol.modificationStates.append(mod)
	self.modificationSites.append( (modificationName, modificationArray[0]))
	
    def parseMolLine(self, theLine):
    
	moleculeString = theLine.strip().split()[-1]
	
	nameIndex = moleculeString.find('(')
	if nameIndex == -1:
	    self.setName(moleculeString.strip())
	    return

	moleculeName, bindings = moleculeString[:nameIndex], moleculeString[nameIndex:]
	self.setName(moleculeName)
	bindings = bindings[:-1]
	bindings = bindings[1:]


	prebindingsArray = bindings.split(',')
	
	bindingSites = [x for x in prebindingsArray if "~" not in x]
	modificationSites = [x for x in prebindingsArray if "~" in x]
	for x in bindingSites:
	    self.addBindingSite(x)

	for x in modificationSites:
	    self.addModificationSite(x)

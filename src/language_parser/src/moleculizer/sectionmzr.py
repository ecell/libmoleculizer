###############################################################################
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    
#
###############################################################################

from parsedline import ParsedLine
from xmlobject import XmlObject
import pdb

class MoleculizerSection( object ) :
    translation_mode = "STRICT"

    def __init__(self, sectionName, sectionBlock):
        self.section_name = sectionName
        self.original_block = sectionBlock[:]
        self.parsed_lines = [ ParsedLine(line) for line in self.original_block ]

    def getMolName(self, molName, ndx):
        return molName + "-" + str(ndx)
        
    def getParsedLines(self):
        return self.parsed_lines

    def getSectionName(self):
        return self.section_name

    def getOriginalBlock(self):
        return self.original_block[:]
        
    def writeParsedComplexAsPlex( self, complex, xmlobject):
        plexElmt = XmlObject( "plex", xmlobject)

        for ndx in range(len(complex.getMols())):
            parsedMol = complex.getMols()[ndx]

            molInst = XmlObject("mol-instance", plexElmt)
            molInst.addAttribute( "name", self.getMolName( parsedMol.getName(), ndx ))
            
            molTypeElmt = XmlObject("mol-ref", molInst)
            molTypeElmt.addAttribute( "name", parsedMol.getName() )

        for bindingID in complex.getBindings().keys():
            # Find the two bindings with that 
            newBindingElmt = XmlObject("binding", plexElmt)

            try:
                molNdx1, molNdx2 = complex.getBindings()[bindingID]
            except ValueError:
                pdb.set_trace()
                raise Exception("Error unpacking bindings '%s' associated with complex '%s' -- Binding ID only has one value." % (repr(complex.getBindings()[bindingID]), complex.getOriginalLine(), bindingID))
            
            parsedMol1 = complex.getMols()[molNdx1]
            parsedMol2 = complex.getMols()[molNdx2]
            bindingSite1 = ""
            bindingSite2 = ""

            if parsedMol1.isModMol():
                for bindingSite in parsedMol1.getBindingSiteList():
                    parsedBindingSite = parsedMol1.getBindingSiteWithName( bindingSite )
                    if parsedBindingSite.hasBindingToken(bindingID):
                        bindingSite1 = parsedBindingSite.getName()
                        break
                else:
                    print "Error, binding site 1 not found"
                    raise Exception()
            else:
                bindingSite1 = parsedMol1.getName()
         
            if parsedMol2.isModMol():
                for bindingSite in parsedMol2.getBindingSiteList():
                    bindingSite = parsedMol2.getBindingSiteWithName( bindingSite )
                    if bindingSite.hasBindingToken( bindingID ):
                        bindingSite2 = bindingSite.getName()
                        break
                else:
                    print "Error, binding site 2 not found"
                    raise Exception()
            else:
                bindingSite2 = parsedMol2.getName()
            
            molInstance1 = XmlObject("mol-instance-ref", newBindingElmt)
            molInstance1.addAttribute( "name", self.getMolName(parsedMol1.getName(), molNdx1))
            
            bindingSiteRef1 = XmlObject("binding-site-ref", molInstance1)
            bindingSiteRef1.addAttribute("name", bindingSite1)


            molInstance2 = XmlObject("mol-instance-ref", newBindingElmt)
            molInstance2.addAttribute( "name", self.getMolName(parsedMol2.getName(), molNdx2))
            
            bindingSiteRef2 = XmlObject("binding-site-ref", molInstance2)
            bindingSiteRef2.addAttribute("name", bindingSite2)

        return plexElmt
        

#     def __molHasBindingToken( self, parsedMol, bindingToken):
#         for parsedBnd in parsedMol.getBindingSiteList():
#             if parsedBnd.hasBindingToken( bindingToken ):
#                 return True
#         else:
#             return False



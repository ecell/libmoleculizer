from sectionmzr import MoleculizerSection
from xmlobject import XmlObject
from section_xcpt import *

class AllosterySection( MoleculizerSection ):
    def __init__(self, name, allosterySection):
        MoleculizerSection.__init__(self, name, allosterySection)
        return

    def getMolName(self, molName, ndx):
        return molName + "-mol-" + str(ndx)

    def assertSanityOfAllosteryLine(self, line):
        if not len(line.getParsedComponents() ) == 1:
            raise Exception()

        if not line.getParsedComponents()[0].isComplex():
            raise Exception()

    def writeParsedAllostericComplexAsPlex( self, complex, xmlobject):
        plexElmt = XmlObject( "plex", xmlobject)

        for ndx in range(len(complex.getMols())):
            parsedMol = complex.getMols()[ndx]

            molInst = XmlObject("mol-instance", plexElmt)
            molInst.addAttribute( "name", self.getMolName( parsedMol.getName(), ndx ))
            
            molTypeElmt = XmlObject("mol-ref", molInst)
            molTypeElmt.addAttribute( "name", parsedMol.getName() )

        for bindingID in complex.getBindings().keys():
            # Find the two bindings with that 
            newBindingElmt = XmlObject("binding-site", plexElmt)
            
            molNdx1, molNdx2 = complex.getBindings()[bindingID]
            
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
            bindingSiteRef2.addAttribute("name", bindingSite1)
            
    def writeAllostericSitesElementToComplex( self, complex, xmlobject):
        return

    def molHasBindingToken( self, parsedMol, bindingToken):
        for parsedBnd in parsedMol.getBindingSiteList():
            if parsedBnd.hasBindingToken( bindingToken ):
                return True
        else:
            return False

class AllostericPlexesSection( AllosterySection ):
    def __init__(self, allostericPlexSection):
        AllosterySection.__init__(self, "AllostericPlex", allostericPlexSection)

    def writeAllostericPlexesSection(self, xmlobject):
        for line in self.getParsedLines():
            self.assertSanityOfAllosteryLine(line)
            self.writeAllostericPlexLineToXml( line, xmlobject)

    def writeAllostericPlexLineToXml( self, lineToWrite, xmlobject):
        alloPlexElmt = XmlObject("allosteric-plex", xmlobject)

        self.writeParsedAllostericComplexAsPlex( lineToWrite.getParsedComponents()[0], alloPlexElmt)
        # self.writeAllostericSitesElementToComplex( lineToWrite.getParsedComponents()[0], alloPlexElmt)

        return
        
class AllostericOmnisSection( AllosterySection ):
    def __init__(self, allostericOmnisSection):
        AllosterySection.__init__(self, "AllostericOmni", allostericOmnisSection)

    def writeAllostericOmnisSection(self, xmlobject):
        for line in self.getParsedLines():
            self.assertSanityOfAllosteryLine(line)
            self.writeAllostericOmniLineToXml( line, xmlobject)

    def writeAllostericOmniLineToXml( self, lineToWrite, xmlobject):
        alloOmniElmt = XmlObject("allosteric-omni", xmlobject)

        self.writeParsedAllostericComplexAsPlex( lineToWrite.getParsedComponents()[0], alloOmniElmt )

        # self.writeAllostericSitesElementToComplex( lineToWrite.getParsedComponents()[0], alloOmniElmt )
        
        return

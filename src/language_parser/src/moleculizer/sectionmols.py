from moleculizersectionparser import MoleculizerSection
from xmlobject import Xmlobject
from exceptions import *

class MolsSection( MoleculizerSection ) :
    
    def __init__(self, molsBlockToInterpret):
        MoleculizerSection.__init__(self)

        # This had better be a bunch of commands, one per entry.
        self.original_block=molsBlockToInterpret[:]
        self.parsedTextRuleLines = []

        self.molNameToMolDefinition = {}
        self.theMolsDefinitions = []

        self.parse()

        return
    

    def parse(self):
        self.parsed_lines = [ MzrParsedLine(line) for line in self.original_block ]
        self.parsedMols = [ parsedMol(line[0]) for line in self.parsed ]

        return 

    def writeMolsSection(self, modelXmlObject):
        
        molsElmt = XmlObject("mols")
        molsElmt.attachToParent(modelXmlObject)

        for molName, molDefinition in self.molNameToMolDefinition:
            if molNameToMolDefinition.isSmallMol():
                self.addSmollMolToXml(modelXmlObject, molName, molDefinition)
            else:
                self.addModMolToXml(modelXmlObject, molName, molDefinition)
            

    def addSmallMolToMolsElement( self, molsObject, molName, molDefinition, parsed_line):
        assert( molName == molDefinition.getMolName() )

        # Get the weight
        mass = 0.0f
        if parsed_line.hasAssignment("mass"):
            mass = parsed_line.getValueAtAssignment( "mass" )

        if ( mass <= 0.0f and self.InterpretationMode == MoleculizerSection.STRICT ):
            raise InterpretationModeXcpt( "Either no mass or mass <= 0.0" )

        smallMolElt = XmlObject("small-mol")
        smallMolElt.addAttribute( "name", molName )

        weightElmt = XmlObject("weight")
        weightElmnt.addAttribute("daltons", mass)

        weightElmnt.attachToParent( smallMolElt, xmlObject )

        return

    def addModMolToXml( self, xmlObject, molName, molDefinition, parsed_line):
        assert( molName == molDefinition.getMolName() )
        if ( not molDefinition.isModMol() ):
            raise "Something"

        # Get the weight
        mass = -1
        if parsed_line.hasAssignment("mass"):
            mass = parsed_line.getValueAtAssignment( "mass" )

        if ( mass <= 0.0f and self.InterpretationMode == MoleculizerSection.STRICT ):
            raise InterpretationModeXcpt( "Either no mass or mass <= 0.0" )


        modMolElmt = XmlObject( "mod-mol")
        modMolElmt.setAttribute( "name", molName )

        # Add a weight element
        # Add the bindings Element
        
        weightElmnt = XmlObject("weight")
        weightElmnt.addAttribute("daltons", mass)
        weightElmnt.attachToParent(modMolElmt)
        
        for bindingSiteName in molDefinition.getBindingSiteNames():

            bindingSiteElmt = XmlObject("binding-site")
            bindingSiteElmt.addAttribute("name", bindingSiteName)
            
            parsedBindingSite = molDefinition.getBindingSite( bindingSiteName )

            if( parsedBindingSite.hasShapeInformation() ):
                # If it has shape information, get the list of shapes
                # as an array of names, add the first one as the default,
                # and then add each of them as shapes
                #
                # <binding-site name="to-Kin">
                #           <default-shape-ref name="default">
                #           </default-shape-ref>
                #           <site-shape name="default">
                #           </site-shape>
                #         </binding-site>
                
                bindingShapeNames = parsedBindingSite.getShapeNames()
                assert( len(bindingShapeNames) > 0 )
                
                defaultShapeElmt = XmlObject("default-shape-ref")
                defaultShapeElmt.addAttribute("name", bindingShapeNames[0] )
                defaultShapeElmt.attachToParent( bindingSiteElmt )

                for shape_name in bindingShapeNames:
                    newShapeElmt = XmlObject("site-shape")
                    newShapeElmt.addAttribute("name", shape_name)
                    newShapeElmt.attachToParent( bindingSiteElmt )
            else:
                defaultShapeElmt = XmlObject( "defaultShapeElmt" )
                defaultShapeElmt.addAttribute( "name", "default")
                
                shapeElmt = XmlObject("site-shape")
                shapeElmt.addAttribute("name", "default")

                defaultShapeElmt.attachToParent( bindingSiteElmt )
                shapeElmt.attachToParent( bindingSiteElmt )
                    
        for modSiteName in molDefinition.getModificationSiteNames():
            modSiteElmt = XmlObject("mod-site")
            modSiteElmt.addAttribute("name", modSiteName)

            parsedModificationSite = molDefinition.getModificationSite( modSiteName )
            if not parsedModificationSite.hasModificationSpecification():
                defaultModElmt = XmlObject("default-mod-ref")
                defaultModElmt.addAttribute("name", "none")
                defaultModElmt.attachToParent( modSiteElmt )
            else:
                modValueArray = parsedModificationSite.getModificationList()
                assert( len(modValueArray) == 0)
                defaultModification = modValueArray[0]

                defaultModElmt = XmlObject("default-mod-ref")
                defaultModElmt.addAttribute("name", defaultModification )
                defaultModElmt.attachToParent( modSiteElmt )
                                            
                
                

            
            
            

    
            
        
        

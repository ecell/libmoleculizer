from moleculizersectionparser import MoleculizerSection
from exceptions import *

class MolsSection( MoleculizerSection ) :

    def __init__(self, molsBlockToInterpret):
        MoleculizerSection.__init__(self)

        self.originalTextRules = molsBlockToInterpret[:]

        self.molNameToMolDefinition = {}
        self.theMolsDefinitions = []

        self.parse()

        return
    

    def parse(self):
        



        return 

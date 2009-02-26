from moleculizermol import MoleculizerMol

class ParsedComplex( ParsedMzrMoleculeToken ):
    def __init__(self, complexTxt):
        ParsedMzrMoleculeToken.__init__(complexTxt)
        
        parsedMolsInComplex = [] 
        parsedBindingsInComplex = []

    def getNumberMolsInComplex(self):
        return len(parsedMolsInComplex)

    def getNumberBindingsInComplex(self): 
        return len(parsedBindingsInComplex)

    def getMols(self):
        return parsedMolsInComplex

    def getBindings(self):
        return parsedBindingsInComplex

     def performSetOfSanityChecks(self):
         pass

     parsed_mols = property( getMols )
     parsed_bindings = property( getBindings )
    

        

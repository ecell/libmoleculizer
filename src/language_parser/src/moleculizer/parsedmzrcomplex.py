from parsedmzrmol import MoleculizerMol

class ParsedComplex( ParsedMzrMoleculeToken ):
    def __init__(self, complexDefinition):
        ParsedMzrMoleculeToken.__init__(complexDefinition)
        
        parsedMolsInComplex = [] 
        parsedBindingsInComplex = []

        self.__parse()

    def __parse(self):
        # For now, the '.' char seperated mols
        self.parsedMolsInComplex = [ ParsedMoleculizerMol(x) for x in self.line.split(".")]
        self.parsedBindingsInComplex = createArrayOfParsedBindings(x)
        
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
    

        

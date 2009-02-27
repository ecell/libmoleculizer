import pdb
from parsedmzrtoken import ParsedMzrToken
from parsedmzrmol import ParsedMoleculizerMol

class ParsedComplex( ParsedMzrToken ):
    def __init__(self, complexTxt):
        ParsedMzrToken.__init__(self, complexTxt)

        self.parsed_mols_in_complex = [] 
        self.parsed_bindings_in_complex = {}

        self.parse()

    def parse(self):
        complexToken = self.getOriginalLine()
        individualMolTokens = complexToken.split(".")

        for molToken in individualMolTokens:
            self.parsed_mols_in_complex.append(ParsedMoleculizerMol(molToken))

        for parsedMolNdx in range(len(self.parsed_mols_in_complex)):
            parsedMol = self.parsed_mols_in_complex[parsedMolNdx]
            if parsedMol.representsBoundMolecule():
                # Register it's binding sites
                if parsedMol.isSmallMol():
                    # It will only have one
                    bindingID = parsedMol.getSmallMolBindingToken()
                    self.__addMolNdxToBoundMap( bindingID, parsedMolNdx)
                elif parsedMol.isModMol():
                    for bindingSiteName in parsedMol.getBindingSiteList():

                        bindingSite = parsedMol.getBindingSiteWithName( bindingSiteName )

                        if bindingSite.hasBindingSiteSpecification() and bindingSite.getBindingSiteSpecification().hasHalfBindingSpecification():
                            bindingID = bindingSite.getBindingSiteSpecification().getHalfBindingSpecification().getBindingToken()
                            self.__addMolNdxToBoundMap( bindingID, parsedMolNdx)
                else:
                    raise Exception()


        pdb.set_trace()
        
        connectivityCheck = [ len(self.parsed_bindings_in_complex[k]) for k in self.parsed_bindings_in_complex.keys() ]

        assert( len(connectivityCheck)== 0 or (min(connectivityCheck) == max(connectivityCheck) == 2) )

    def getNumberMolsInComplex(self):
        return len(self.parsed_mols_in_complex)

    def getNumberBindingsInComplex(self): 
        return len(self.parsed_bindings_in_complex)

    def getMols(self):
        return self.parsed_mols_in_complex

    def getBindings(self):
        return self.parsed_bindings_in_complex

    def __addMolNdxToBoundMap(self, bindingNdx, molNdx):
        if bindingNdx not in self.parsed_bindings_in_complex.keys():
            self.parsed_bindings_in_complex[ bindingNdx ] = []

        self.parsed_bindings_in_complex[bindingNdx].append( molNdx )


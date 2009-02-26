from parsedbindingsite.py import ParsedBindindSite
from parsedmodificationsite.py import ParsedBindindSite

class ParsedMoleculizerMol( ParsedMzrToken ):
    def __init__(self, molDefinition ):
        ParsedMzrToken.__init__(self, molDefinition )

        # I want to be strict about forbidding bindings
        # Parsing of complexes must factor them out first

        self.isSmallMol = True

        self.molName = ""

        self.theParsedBindings = []
        self.theParsedModificationSites = []

        self.parse()

    def isSmallMol(self):
        return self.isSmallMol;

    def isModMol(self):
        return not self.isSmallMol

    def getBindingList(self):
        if self.isSmallMol():
            raise Exception()

        return [parsed_binding.getName() for parsed_binding in self.theParsedBindings]

    def getModificationList(self):
        if self.isSmallMol():
            raise Exception()

        return [ parsed_modification.getName() for parsed_modification in self.theParsedModificationSites]
        
    def getBindingSiteWithName(self, bindingName ):
        # Should I return the one bindng
        if self.isSmallMol():
            raise Exception()

        return

    def getModificationSiteWithName(self, modificationSite ):
        if self.isSmallMol():
            raise Exception()

        return 

    def parse(self):

        # Token here means a single mol definition alpha, or alpha(first,second), or whatever.
        tokenToParse = self.getOriginalLine()
        
        if not self.pre_parse_sanity_checks():
            raise ParsedMolSyntacticalException( tokenToParse )

        if "(" in tokenToParse:
            self.__parseModMol()
        else:
            self.__parseSmallMol()

        if not self.post_parse_sanity_checks:
            raise ParsedMolSemanticException( self.getOriginalLine() )


    def __parseSmallMol(self):
        # Parse 
        smallMolTokenToParse = self.getOriginalLine()
        
        self.molName = smallMolTokenToParse

    def __parseModMol(self):
        # The Token is a Mod Mol, parse it
        
        modMolTokenToParse = self.getOriginalLine()
        assert( modMolTokenToParse.count("(") == 1 and modMolTokenToParse.count(")") ==2 )

        toks = modMolTokenToParse.split("(")
        self.molName = toks[0]
        
        structTokensStr = toks[1][:-1] # The last charector is the last ")" at the end of the mol.
        
        structTmpTokens = structTokensStr.split(",")
        structTokens = []
        
        for struct_token in structTokens:
            if struct_token.starts_with("*"):
                self.theParsedModificationSites.append( ParsedModificationSite(struct_token) )
            else:
                self.theParsedBindings.append( ParsedBindingsSite( struct_token ))
        


        
        
        

class ParsedModificationSite( ParsedMzrToken ):
    def __init__(self, modificationSiteToken):
        ParsedMzrToken.__init__(self, modificationSiteToken)
        
        assert(modificationSiteToken.startswith("*"))

        self.modificationSiteName = ""
        self.definedModificationValueSpecifier = False
        self.modificationSiteSpecificationTokens = []

        self.parse()

        return

    def parse(self):
        return

    def getModificationSiteName(self):
        return self.modificationSiteName
    
    def modificationSiteHasSpecification(self):
        return 
                                             

    


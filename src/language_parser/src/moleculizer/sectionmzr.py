###############################################################################
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    
#
###############################################################################

from parsedline import ParsedLine

class MoleculizerSection( object ) :
    translation_mode = "STRICT"

    def __init__(self, sectionName, sectionBlock):
        self.section_name = sectionName
        self.original_block = sectionBlock[:]
        self.parsed_lines = [ ParsedLine(line) for line in self.original_block ]
        
    def getParsedLines(self):
        return self.parsed_lines

    def getSectionName(self):
        return self.section_name

    def getOriginalBlock(self):
        return self.original_lines[:]
        
        



class ParsedMzrMoleculeToken( object ):
    def __init__(self, original_line, line_number =- 1):
        self.original_line = original_line
        self.line_number = line_number

    def getLineNumber(self):
        return line_number

    line_number = property( getLineNumber )
            
        

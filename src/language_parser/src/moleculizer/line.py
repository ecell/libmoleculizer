import pdb

class MzrRuleFileLine( object):
    def __init__(self, lineTxt):
        # A MzrRuleFileLine is a list of tokens (each of which can be a molDefinition,
        # 
        self.lineText = lineTxt
        self.tokens = []

        self.__parse()

    def __parse(self):
        tmpLine = self.lineText

        tmp_unprocessed_tokens = tmpLine.spit(",")
        tmp_unprocessed_tokens = [x.strip() for x in tmp_unprocessed_tokens]

        self.tokens = getBalancedArray( tmp_unprocessed_tokens )


def getBalancedArray( arrayToBalance ):
    if not EachEntryIsParenBalanced( arrayToBalance ):
    # Combine the ..., ndx_i, ndx_(i+1) where ndx_i is the smallest i not balanced
        return getBalancedArray( GetIncrementallyBetterArray( arrayToBalance ) )
    else:
        return arrayToBalance

def GetIncrementallyBetterArray( anArray ):

    values = [ StringIsParenBalenced(x) for x in anArray]

    assert( False in values )

    badNdx = values.index( False )
    
    combinedTokens = anArray[badNdx] + anArray[badNdx + 1]

    returnArray = anArray[ : badNdx]
    returnArray.append( combinedTokens )
    returnArray.extend( anArray[badNdx + 2 : ] )
    
    return returnArray
    
def EachEntryIsParenBalanced( array ):
    entries = [ StringIsParenBalenced(x) for x in array ]
    
    returnVal = True
    for val in entries:
        returnVal = returnVal and val

    return returnVal

def StringIsParenBalenced(line):
    return ( line.count("(") == line.count(")") and 
             line.count("[") == line.count("]") and 
             line.count("{") == line.count("}") )

        



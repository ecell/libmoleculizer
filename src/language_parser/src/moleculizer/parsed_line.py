###############################################################################
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    
#
###############################################################################

import pdb

from util import *

class MzrRuleFileLine( object):
    def __init__(self, lineTxt):
        # A MzrRuleFileLine is a list of tokens (each of which can be a molDefinition,
        # 
        self.lineText = lineTxt
        self.tokens = []

        self.__parse()

    def __parse(self):
        # We can assume that self.line text is a single comment, with a semi-colon ending it
        
        # No semi-colon
        self.tmpLine = self.snip_statement(self.lineText)
        
        tmpLine = self.lineText

         
        tmp_unprocessed_tokens = tmpLine.spit(",")
        tmp_unprocessed_tokens = [x.strip() for x in tmp_unprocessed_tokens]
         
        self.tokens = getBalancedArray( tmp_unprocessed_tokens )



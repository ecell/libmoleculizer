###############################################################################
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    
#
###############################################################################


from sectionmzr import MoleculizerSection
from xmlobject import XmlObject
from section_xcpt import *

class DimerGenSection( MoleculizerSection ):
    def __init__(self, dimerBlock):
        MoleculizerSection.__init__(self, "DimerGenSection", dimerBlock)
        return 

class OmniGenSection( MoleculizerSection ):
    def __init__(self, omniBlock):
        MoleculizerSection.__init__(self, "OmniGenSection", omniBlock)
        return 

class UniGenSection( MoleculizerSection ):
    def __init__(self, uniBlock):
        MoleculizerSection.__init__(self, "UniGenSection", uniBlock)
        return 


class ReactionRulesSection( MoleculizerSection ):
    def __init__(self, rxnBlock, dimerGenBlock, omniGenBlock, uniGenBlock):
        MoleculizerSection.__init__(self, "ReactionRulesSection", rxnBlock)

        self.the_dimer_section = DimerGenSection( dimerGenBlock )
        self.the_omni_gen_section = OmniGenSection( omniGenBlock )
        self.the_uni_gen_secont = UniGenSection( uniGenBlock )
        return



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

class ExplicitSpeciesSection( MoleculizerSection ):
    def __init__(self, explicitSpeciesBlock ):
        MoleculizerSection.__init__(self, "ExplicitSpeciesSection", explicitSpeciesBlock)
        return
        

###############################################################################
# BNGMZRConverter - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
#
# Moleculizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Moleculizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Moleculizer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Original Author:
#   Nathan Addy, Scientific Programmer	Email: addy@molsci.org
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#                     
#   
###############################################################################

def getMolNameFromFullMolSpec( bngMolSpecification ):
    # This function takes a MOL name - Ste11(bnd_1, bnd_2, a~0~1) and returns
    # the name, in this case Ste11
    
    bngMolSpecification = bngMolSpecification.strip()

    ndx = bngMolSpecification.find('(')
    if ndx == -1:
        return bngMolSpecification
    else:
        return bngMolSpecification[:ndx]


def getBindingSpecificationFromNdx(bngComplexSpecification, bindingNdx):
    # This function takes a complex name like Ste5(ste_11!1, ste_12, ste13).Ste11(ste_5!1, b~0~1)
    # and an int like 1 and reurns [(Ste5, ste11), (Ste11, ste_5)].
    
    bindingNdxStr = str(bindingNdxStr)
    
    bngMolComponentList = bngComplexSpecification.split('.')
    bindingComponents = []
    
    for component in bngMolecularComponentList:
        if bindingNdxStr in component:
            molName = getMolNameFromMolSpec( component )
            bindingName = getBindingNameWithBoundID( component, bindingNdxStr )
            bindingComponents.append( (molName, bindingName) )
    else:
        assert( len(bindingComponents) == 2)
        bindingComponents.sort()
        return bindingComponents
    raise ValueError, "Error, bad binding.  Plus I need a better error message..."

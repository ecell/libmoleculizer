###############################################################################
# BNGMZRConverter - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007, 2008 The Molecular Sciences Institute
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Contact information:
#   Nathan Addy, Scientific Programmer	Voice: 510-981-8748
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#   2168 Shattuck Ave.                  
#   Berkeley, CA 94704
###############################################################################

import util
import pdb

class MatchObject:
    # This is an object that matches another.  It is the object upon which all reactions
    # are sorted.  Any two reactions which "sort" the same are allosteric versions of
    # identical reactions.

    # So we are sorting lexicographically on (self.dimerization, bindingAffected), is
    # basically what this thing does.

    # Adding in modification reactions...
    def __init__(self, reactants, products):
        self.numberOfReactants = len(reactants)

        self.reactantList = reactants[:]
        self.productList = products[:]

        self.modificationReaction = 0
        self.dimerization = 0
        self.decomposition = 0
        self.significantBinding = []

        try:
            self.initialize()
        except:
            # Don't quite know what to do here.  Can we raise exceptions properly in a constructor?
            pass

    def initialize(self):

        if len(self.reactantList) == 1 and len(self.productList) == 2:
            self.decomposition, self.dimerization, self.modificationReaction = True, False, False
        elif len(self.reactantList) == 2 and len(self.productList) == 1:
            self.decomposition, self.dimerization, self.modificationReaction = False, True, False
        elif len(self.reactantList) == 1 and len(self.productList) == 1:
            self.decomposition, self.dimerization, self.modificationReaction = False, False, True
        else:
            raise ValueError, "Reaction not properly formed in MatchObject initializer..."

        if self.modificationReaction:
            self.significantBinding = util.calculateSignificantModificationSite(self.reactantList, self.productList)
        else:
            try:
                self.significantBinding = util.calculateSignificantBinding(self.reactantList, self.productList)
            except:
                # Whoops bad reaction thingy.
                print "Poorly formed reaction..."

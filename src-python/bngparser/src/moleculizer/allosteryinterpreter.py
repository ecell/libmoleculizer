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

from namegenerator import NameGenerator
from xmlobject import XmlObject
import util
import pdb

class AllosteryInterpreter:
    def __init__(self, moleculeInterpreter):
        self.moleculeInterpreter = moleculeInterpreter
        self.allosteryMapping = {}
        self.bindingAllosteryNameGenerator = {}

    def addAllostericOmnis(self, allostericOmnisElement):
        # pdb.set_trace()
        for bindingWithAllostery in self.allosteryMapping.keys():
            owningMolName, allostericBindingSiteName = bindingWithAllostery
            complexAllosteryPairs = [(key, self.allosteryMapping[bindingWithAllostery][key]) for key in self.allosteryMapping[bindingWithAllostery].keys()]
            for complex_allostery_pair in complexAllosteryPairs:
                complex, allostericName = complex_allostery_pair
                if allostericName == "default":
                    continue

                allostericOmniElmt = XmlObject("allosteric-omni")
                allostericOmnisElement.addSubElement(allostericOmniElmt)

                # Now prepare the four sub-elements
                plexElmt = XmlObject("plex")
                instanceStatesElmt = XmlObject("instance-states")
                unboundSitesElmt = XmlObject("unbound-sites")
                boundSitesElmt = XmlObject("bound-sites")
                allostericStatesElmt = XmlObject("allosteric-sites")

                # First we construct the molRefToInstanceDict
                molRefToInstanceDict = {}
                tmpMolInstanceInt = {}

                molArray = complex.split('.')

                for molNdx in xrange(len(molArray)):
                    mol = molArray[molNdx]
                    parenNdx = mol.find("(")
                    molType = mol[:parenNdx]

                    if molType not in tmpMolInstanceInt:
                        tmpMolInstanceInt[molType] = 1
                        molInstance = "the-" + molType
                    else:
                        molInstance = molType + "-" + str(tmpMolInstanceInt[molType])
                        tmpMolInstanceInt[molType] = tmpMolInstanceInt[molType] + 1

                    molRefToInstanceDict[molInstance] = molType
                    molArray[molNdx] = molInstance + mol[parenNdx:]

                # So now we calculate -- the plex.
                bindingDict = {}
                genericBoundSitesArray = []
                unboundSitesArray = []
                modificationSitesArray = []
                for molSpec in molArray:
                    # Process the molNames
                    molInstance = util.getMolNameFromFullMolSpec(molSpec)
                    molInstanceElmt = XmlObject("mol-instance")
                    molInstanceElmt.addAttribute("name", molInstance)
                    molRefElmt = XmlObject("mol-ref")
                    molRefElmt.addAttribute("name", molRefToInstanceDict[molInstance])
                    molRefElmt.attachToParent(molInstanceElmt).attachToParent(plexElmt)

                    # Collect the bindings
                    bindingSites = util.getBindingSites(molSpec)

                    explicitlyBoundBindingSites = [x.split("!") for x in bindingSites if ("!" in x and x[-1] != "!" and x[-1] != "+")]
                    boundBindingSites = [x.split("!") for x in bindingSites if ("!" in x and x[-1] != "!" and x[-1] == "+")]
                    unboundBindingSites = [x for x in bindingSites if "!" not in x]

                    modificationSites = [x.split("~") for x in util.getModificationList(molSpec)]

                    for halfBinding in explicitlyBoundBindingSites:
                        bindingSite, bindingNdx = halfBinding[0], int(halfBinding[1])
                        if not bindingNdx in bindingDict:
                            bindingDict[bindingNdx] = []
                        bindingDict[bindingNdx].append( (molInstance, bindingSite))

                    # Collect the unbound binding sites.
                    for unboundBindingSite in unboundBindingSites:
                        unboundSitesArray.append( (molInstance, unboundBindingSite) )

                    for boundSite in boundBindingSites:
                        genericBoundSitesArray.append( (molInstance, boundSite[0]))

                    for modificationSite in modificationSites:
                        modificationSitesArray.append( (molInstance, modificationSite[0], modificationSite[1]))

                # Do the bindings
                listOfBindings = [ bindingDict[key] for key in bindingDict.keys()]

                for binding in listOfBindings:
                    try:
                        assert(len(binding)==2)
                    except:
                        pdb.set_trace()
                        a = 1
                for binding in listOfBindings:
                    bindingElmt = XmlObject("binding")
                    bindingElmt.attachToParent(plexElmt)
                    for halfBinding in binding:
                        molInstance, bindingRef = halfBinding
                        molInstanceRefElmt = XmlObject("mol-instance-ref")
                        molInstanceRefElmt.addAttribute("name", molInstance)
                        bindingSiteRefElmt = XmlObject("binding-site-ref")
                        bindingSiteRefElmt.addAttribute("name", bindingRef)
                        bindingSiteRefElmt.attachToParent(molInstanceRefElmt).attachToParent(bindingElmt)

                # Now we add in the modificationSitesArray
                for modMolName in set( [ x[0] for x in modificationSitesArray]):
                    modMolModifications = [x for x in modificationSitesArray if x[0] == modMolName]
                    modMolInstanceRefElmt = XmlObject("mod-mol-instance-ref")
                    modMolInstanceRefElmt.attachToParent(instanceStatesElmt)
                    modMolInstanceRefElmt.addAttribute("name", modMolName)
                    modMapElmt = XmlObject("mod-map")
                    modMapElmt.attachToParent(modMolInstanceRefElmt)
                    for modification in modMolModifications:
                        modSiteRefElmt = XmlObject("mod-site-ref")
                        modSiteRefElmt.addAttribute("name", modification[1])
                        modRefElmt = XmlObject("mod-ref")
                        modRefElmt.addAttribute("name", modification[2])
                        modRefElmt.attachToParent(modSiteRefElmt).attachToParent(modMapElmt)

                # Now we add in the unboundSitesArray
                for molWithUnboundBindings in set( [ x[0] for x in unboundSitesArray]):
                    unboundSites = [x for x in unboundSitesArray if x[0] == molWithUnboundBindings]
                    for site in unboundSites:
                        instanceRefElmt = XmlObject("instance-ref")
                        instanceRefElmt.addAttribute("name", molWithUnboundBindings)
                        siteRefElmt = XmlObject("site-ref")
                        siteRefElmt.addAttribute("name", site[1])
                        siteRefElmt.attachToParent(instanceRefElmt).attachToParent(unboundSitesElmt)

                # Now we add in the boundSitesArray
                for molWithGenericBoundBindings in set([ x[0] for x in genericBoundSitesArray]):
                    genericBoundSites = [x for x in genericBoundSitesArray if x[0] == molWithGenericBoundBindings]
                    for site in genericBoundSites:
                        instanceRefElmt = XmlObject("instance-ref")
                        instanceRefElmt.addAttribute("name", molWithGenericBoundBindings)
                        siteRefElmt = XmlObject("site-ref")
                        siteRefElmt.addAttribute("name", site[1])
                        siteRefElmt.attachToParent(instanceRefElmt).attachToParent(boundSitesElmt)
                        


                molInstanceRefElmt = XmlObject("mol-instance-ref")
                molInstanceRefElmt.addAttribute("name", "the-" + owningMolName)
                bindingSiteRefElmt = XmlObject("binding-site-ref")
                bindingSiteRefElmt.addAttribute("name", allostericBindingSiteName)
                siteShapeRefElmt = XmlObject("site-shape-ref")
                siteShapeRefElmt.addAttribute("name", allostericName)
                siteShapeRefElmt.attachToParent(bindingSiteRefElmt).attachToParent(molInstanceRefElmt).attachToParent(allostericStatesElmt)

                ## End of function
                plexElmt.attachToParent(allostericOmniElmt)
                if instanceStatesElmt.getNumberOfChildren() > 0:
                    instanceStatesElmt.attachToParent(allostericOmniElmt)
                if unboundSitesElmt.getNumberOfChildren() > 0:
                    unboundSitesElmt.attachToParent(allostericOmniElmt)
                if boundSitesElmt.getNumberOfChildren():
                    boundSitesElmt.attachToParent(allostericOmniElmt)
                allostericStatesElmt.attachToParent(allostericOmniElmt)
                
                
    def DEBUG_getAllosteryMappingSize(self):
        self.diff = []
        try:
            self.runOnce
        except:
            # This is the first time this function has been run
            self.lastChecked = []
            self.runOnce = 0

        size = 0
        for key in self.allosteryMapping.keys():
            size += len(self.allosteryMapping[key])

        for key in self.allosteryMapping.keys():
            for entry in self.allosteryMapping[key]:
                if not (key, entry) in self.lastChecked:
                    self.diff.append( (key, entry))

        self.lastChecked.extend(self.diff)
        return size, self.diff[:]
                      

    def createNewAllosteryAtSite(self, site, complex):
        if not site in self.bindingAllosteryNameGenerator.keys():
            self.bindingAllosteryNameGenerator[site] = NameGenerator(site[0] + "_" + site[1])

        if site not in self.allosteryMapping.keys():
           self.allosteryMapping[site] = {}

        if complex in self.allosteryMapping[site].keys():
            # This complex has already been recorded as an allosteric state of site.
            return 

        if self.checkComplexIsDefault(site, complex):
            allosteryName = "default"
        else:
            allosteryName = self.bindingAllosteryNameGenerator[site].generateNewName()



        self.allosteryMapping[site][complex] = allosteryName
        
        
    def checkComplexIsDefault(self, binding, complex):
        try:
            molName, bindingSiteName = binding
        except:
            pdb.set_trace()
            a = 1
        if "." in complex:
            return False
        if not complex.startswith(molName):
            pdb.set_trace()
            raise ValueError, "Bad complex"
        if "!" in complex:
            return False
        theBindingSites = util.getBindingSites(complex)
        if len(theBindingSites) != 1:
            return False
        for modSite in util.getModificationList(complex):
            modSite = modSite.split("~")
            if modSite[1] != "none":
                return False
        else:
            return True
            

        

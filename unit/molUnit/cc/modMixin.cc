/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001  Walter Lawrence (Larry) Lok.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#include "mol/modMixin.hh"

namespace bnd
{
  double
  modStateMixin::totalWeightDelta(void) const
  {
    double totalDelta = 0.0;
    
    int modNdx = size();
    while(0 < modNdx--) totalDelta += (*this)[modNdx]->weightDelta;

    return totalDelta;
  }

  // Gives the number of modification sites with the given modification.
  int
  modStateMixin::modCount(const modification* pModToCount) const
  {
    int modCount = 0;

    int modNdx = size();
    while(0 < modNdx--)
      {
	if((*this)[modNdx] == pModToCount) ++modCount;
      }

    return modCount;
  }

  // Same as above, but only counts sites that are masked true.
  int
  modStateMixin::modCount(const modification* pModToCount,
			  const std::vector<bool>& rMask) const
  {
    int modCount = 0;

    int modNdx = size();
    while(0 < modNdx--)
      {
	if(rMask[modNdx]
	   && ((*this)[modNdx] == pModToCount)) ++modCount;
      }

    return modCount;
  }

  // Changes a randomly selected site with modification pFromMod to
  // have modification pToMod.  Expect this to be used, e.g. to
  // phosphorylate an unphosphorylated site.
  //
  // Note that this could just as well be used to dephosphorylate
  // a phosphorylated site.
  // 
  // Returns the index of the modified site or -1 if no site had
  // modification fromMod.
  //
  // Do I need this any more???
  // int
  // modStateMixin::pushMod(const modification* pFromMod,
  // 		       const modification* pToMod)
  // {
  //   // Form the inverse cumulative distribution function for the
  //   // uniform distribution on the sites that currently have the
  //   // fromMod modification.
  //   std::vector<int> icdf;
  //   int modNdx = size();
  //   while(0 < modNdx--)
  //     {
  //       if((*this)[modNdx] == pFromMod) icdf.push_back(modNdx);
  //     }

  //   // If any available sites were turned up, select one of them
  //   // at random for modification.
  //   modNdx = -1;
  //   if(0 < icdf.size())
  //     {
  //       modNdx = icdf[samplePop(icdf.size())];
  //       (*this)[modNdx] = pToMod;
  //     }

  //   return modNdx;
  // }

  // Same as above, but acts only on sites masked true.
  // int
  // modStateMixin::pushMod(const modification* pFromMod,
  // 		       const modification* pToMod,
  // 		       std::vector<bool>& rMask)
  // {
  //   std::vector<int> icdf;
  //   int modNdx = size();
  //   while(0 < modNdx--)
  //     {
  //       if(rMask[modNdx]
  // 	 && ((*this)[modNdx] == pFromMod)) icdf.push_back(modNdx);
  //     }

  //   modNdx = -1;
  //   if(0 < icdf.size())
  //     {
  //       modNdx = icdf[samplePop(icdf.size())];
  //       (*this)[modNdx] = pToMod;
  //     }

  //   return modNdx;
  // }

  modMolMixin::modMolMixin
  (const std::map<std::string, const modification*>& rDefaultModMap)
  {
    for(std::map<std::string, const modification*>::const_iterator
	  iModEntry = rDefaultModMap.begin();
	iModEntry != rDefaultModMap.end();
	++iModEntry)
      {
	const std::string& rSiteName = iModEntry->first;
	int siteNdx = modSiteNames.size();

	// Checking uniqueness of modification site names.
	bool insertOk
	  = modSiteNameToNdx.insert(std::make_pair(rSiteName,
						   siteNdx)).second;
	if(insertOk)
	  {
	    modSiteNames.push_back(rSiteName);
	  }
	else
	  {
	    throw duplicateModSiteNameXcpt(rSiteName);
	  }
      }
  }

  class modSubstituter :
    public std::unary_function<std::pair<std::string, const modification*>, void>
  {
    const modMolMixin& rModMol;
    modStateMixin& rTarget;
  public:
    modSubstituter(const modMolMixin& rModMolMixin,
		   modStateMixin& rTargetStateMixin) :
      rModMol(rModMolMixin),
      rTarget(rTargetStateMixin)
    {}

    void
    operator()(const std::pair<std::string, const modification*>& rEntry) const
    {
      const std::string& rSiteName = rEntry.first;
      const modification* pMod = rEntry.second;
      
      std::map<std::string, int>::const_iterator iCatEntry
	= rModMol.modSiteNameToNdx.find(rSiteName);

      // I don't really have a modMol here, just a modMolMixin, so I can't
      // get the name of the mol.
      if(iCatEntry == rModMol.modSiteNameToNdx.end())
	throw unknownModSiteXcpt(rSiteName);

      int modNdx = iCatEntry->second;

      rTarget[modNdx] = pMod;
    }
  };

  modStateMixin modMolMixin::
  substituteModMap(const std::map<std::string, const modification*>& rModMap,
		   const modStateMixin& rSourceStateMixin)
  {
    modStateMixin resultStateMixin(rSourceStateMixin);

    for_each(rModMap.begin(),
	     rModMap.end(),
	     modSubstituter(*this,
			    resultStateMixin));

    return resultStateMixin;
  }

  // modMixin.hh has notes on this.
  modStateMixin modMolMixin::
  indexModMap(const std::map<std::string, const modification*>& rModMap)
  {
    modStateMixin resultStateMixin(0, modSiteNames.size());

    for_each(rModMap.begin(),
	     rModMap.end(),
	     modSubstituter(*this,
			    resultStateMixin));

    return resultStateMixin;
  }

  // modMixin.hh has notes on this.
  bool
  modMolMixin::modStateMatch::operator()(const modStateMixin& rMixinToTest) const
  {
    int modNdx = rMatch.size();
    while(0 < modNdx--) 
      {
	const modification* pMatchMod = rMatch[modNdx];

	if((0 != pMatchMod)
	   && (rMixinToTest[modNdx] != pMatchMod)) return false;
      }
    return true;
  }

  int
  modMolMixin::getModSiteNdx(const std::string& rModSiteName) const
  {
    std::map<std::string, int>::const_iterator iEntry
      = modSiteNameToNdx.find(rModSiteName);

    return iEntry == modSiteNameToNdx.end()
      ? -1
      : iEntry->second;
  }

  // This should be an method of modMol, obviating the stupid pModMol
  // argument.
  int
  modMolMixin::mustGetModSiteNdx(xmlpp::Node* pRequestingNode,
				 const std::string& rModSiteName,
				 modMol* pModMol) const
    throw(unknownModSiteXcpt)
  {
    int ndx = getModSiteNdx(rModSiteName);
    if(ndx < 0) throw unknownModSiteXcpt(pRequestingNode,
					 rModSiteName,
					 pModMol);
    return ndx;
  }
}

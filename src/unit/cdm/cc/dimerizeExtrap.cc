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

#include "cdm/dimerizeExtrap.hh"
#include "cdm/missingDimerizeRateXcpt.hh"
#include "cdm/missingDimerizeInvariantXcpt.hh"

namespace cdm
{
  void
  dimerizeNoExtrap::
  setRate(cpx::siteParam leftParam,
	  cpx::siteParam rightParam,
	  double rate)
  {
    // For the time being, I'm storing the pair in only one order,
    // then searching for both orders when the pair is looked up.
    std::pair<rateMapType::iterator, bool> insertResult
      = rateMap.insert(std::make_pair(std::make_pair(leftParam,
						     rightParam),
				      rate));
    if(! insertResult.second)
      {
	insertResult.first->second = rate;
      }
  }


  double
  dimerizeNoExtrap::
  getRate(const cpx::cxSite<clx::cptPlexSpecies, clx::cptPlexFamily>& rLeftContext,
	  const cpx::cxSite<clx::cptPlexSpecies, clx::cptPlexFamily>& rRightContext) const
  {
    // Now we have to check for both orders in the pair of site shape
    // pointers.
    rateMapType::const_iterator iEntry
      = rateMap.find(std::make_pair(rLeftContext.getSiteParam(),
				    rRightContext.getSiteParam()));
    if(iEntry == rateMap.end())
      {
	iEntry = rateMap.find(std::make_pair(rRightContext.getSiteParam(),
					     rLeftContext.getSiteParam()));
	if(iEntry == rateMap.end())
	  throw missingDimerizeRateXcpt(rLeftContext,
					rRightContext);
      }

    return iEntry->second;
  }

  void
  dimerizeMassExtrap::
  setRate(cpx::siteParam leftParam,
	  cpx::siteParam rightParam,
	  double rate)
  {
    double invariant = fnd::bindingInvariant(rate,
					     leftMass,
					     rightMass);
      
    std::pair<invMapType::iterator, bool> insertResult
      = invariantMap.insert(std::make_pair(std::make_pair(leftParam,
							  rightParam),
					   invariant));
    if(! insertResult.second)
      {
	insertResult.first->second = invariant;
      }
  }

  double
  dimerizeMassExtrap::
  getRate(const cpx::cxSite<clx::cptPlexSpecies, clx::cptPlexFamily>& rLeftContext,
	  const cpx::cxSite<clx::cptPlexSpecies, clx::cptPlexFamily>& rRightContext) const
  {
    // Since the rates were stored with the key pair in only one
    // order, we have to check for both orders when looking it up.
    invMapType::const_iterator iEntry
      = invariantMap.find(std::make_pair(rLeftContext.getSiteParam(),
					 rRightContext.getSiteParam()));
    if(iEntry == invariantMap.end())
      {
	iEntry = invariantMap.find(std::make_pair(rRightContext.getSiteParam(),
						  rLeftContext.getSiteParam()));
	if(iEntry == invariantMap.end())
	  throw missingDimerizeInvariantXcpt(rLeftContext,
					     rRightContext);
      }

    return fnd::bindingRate(iEntry->second,
			    rLeftContext.getPlexWeight(),
			    rRightContext.getPlexWeight());
  }
}

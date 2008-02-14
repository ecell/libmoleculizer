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

#include "dimer/decomposeExtrap.hh"
#include "dimer/missingDecomposeRateXcpt.hh"

namespace dimer
{
  void
  decomposeNoExtrap::
  setRate(cpx::siteParam leftParam,
	  cpx::siteParam rightParam,
	  double rate)
  {
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
  decomposeNoExtrap::
  getRate(const cpx::cxBinding<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rContext) const
    throw(utl::xcpt)
  {
    // Get the shapes (bnd::siteParam) of the sites in the binding.
    std::pair<cpx::siteParam, cpx::siteParam> siteParams
      = rContext.getSiteParams();
      
    // Since the rates were stored with the key pair in only one
    // order, we have to check for both orders when looking it up.
    rateMapType::const_iterator iEntry
      = rateMap.find(siteParams);

    if(iEntry == rateMap.end())
      {
	iEntry = rateMap.find(std::make_pair(siteParams.second,
					     siteParams.first));
	if(iEntry == rateMap.end())
	  throw missingDecomposeRateXcpt(siteParams.first->getName(),
					 siteParams.second->getName());
      }
      
    return iEntry->second;
  }
}

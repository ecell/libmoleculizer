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

#ifndef DECOMPOSEEXTRAP_H
#define DECOMPOSEEXTRAP_H

#include "dimer/dimerXcpt.hh"

namespace dimer
{
  // Base class of rate extrapolators for decomposition reactions.
  class decomposeExtrapolator
  {
  public:
    virtual
    ~decomposeExtrapolator(void)
    {}

    virtual void
    setRate(bnd::siteParam leftParam,
	    bnd::siteParam rightParam,
	    double rate) = 0;
    
    virtual double
    getRate(const plx::bindingInContext& rContext) const = 0;
  };

  // Standard decomposition extrapolator, which doesn't really
  // do anything: the same rate applies regardless of the context.
  class decomposeNoExtrap :
    public decomposeExtrapolator
  {
    typedef
    std::map<std::pair<bnd::siteParam, bnd::siteParam>, double>
    rateMapType;
    
    rateMapType rateMap;

  public:
    // Both for inserting default rates and for writing allosteric
    // rates over default rates.
    // 
    // Storing the rate with the key pair in only one order; this means that
    // search has to look for both orders.
    void
    setRate(bnd::siteParam leftParam,
	    bnd::siteParam rightParam,
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

    // Retrieve decomposition rate for binding in a new complex species.
    double
    getRate(const plx::bindingInContext& rContext) const
    {
      plx::cxBinding cx(rContext);
      
      // Get the shapes (bnd::siteParam) of the sites in the binding.
      std::pair<bnd::siteParam, bnd::siteParam> siteParams
	= cx.getSiteParams();
      
      // Since the rates were stored with the key pair in only one
      // order, we have to check for both orders when looking it up.
      rateMapType::const_iterator iEntry
	= rateMap.find(siteParams);

      if(iEntry == rateMap.end())
	{
	  iEntry = rateMap.find(std::make_pair(siteParams.second,
					       siteParams.first));
	  if(iEntry == rateMap.end())
	    throw missingDecomposeRateXcpt(cx);
	}
      
      return iEntry->second;
    }
  };
    
}

#endif // DECOMPOSEEXTRAP_H

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

#ifndef CDM_DECOMPOSEEXTRAP_H
#define CDM_DECOMPOSEEXTRAP_H

#include "cpx/siteShape.hh"
#include "clx/cptPlexSpecies.hh"
#include "clx/cptPlexFamily.hh"

namespace cdm
{
  // Base class of rate extrapolators for decomposition reactions.
  class decomposeExtrapolator
  {
  public:
    virtual
    ~decomposeExtrapolator(void)
    {}

    virtual void
    setRate(cpx::siteParam leftParam,
	    cpx::siteParam rightParam,
	    double rate) = 0;
    
    virtual double
    getRate(const cpx::cxBinding<clx::cptPlexSpecies, clx::cptPlexFamily>& rContext) const
      throw(utl::xcpt) = 0;
  };

  // Standard decomposition extrapolator, which doesn't really
  // do anything: the same rate applies regardless of the context.
  class decomposeNoExtrap :
    public decomposeExtrapolator
  {
    typedef
    std::map<std::pair<cpx::siteParam, cpx::siteParam>, double>
    rateMapType;
    
    rateMapType rateMap;

  public:
    // Both for inserting default rates and for writing allosteric
    // rates over default rates.
    // 
    // Storing the rate with the key pair in only one order; this means that
    // search has to look for both orders.
    void
    setRate(cpx::siteParam leftParam,
	    cpx::siteParam rightParam,
	    double rate);

    // Retrieve decomposition rate for binding in a new complex species.
    double
    getRate(const cpx::cxBinding<clx::cptPlexSpecies, clx::cptPlexFamily>& rContext) const
      throw(utl::xcpt);
  };
    
}

#endif // CDM_DECOMPOSEEXTRAP_H

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

#ifndef MODKINASEEXTRAP_H
#define MODKINASEEXTRAP_H

#include "plex/cxMolParam.hh"

namespace kinase
{
  // Base class for rate extrapolator for modKinase reactions.
  // This is a binary-type reaction family.
  class modKinaseExtrapolator
  {
  public:
    virtual
    ~modKinaseExtrapolator(void)
    {}

    virtual double
    getRate(const plx::molInContext& rKinaseContext,
	    const plx::molInContext& rSubstrateContext) const = 0;
  };

  // modKinase rate extrapolator that does no extrapolation.  Not
  // recommended; included as an example.
  class modKinaseNoExtrap :
    public modKinaseExtrapolator
  {
    double rate;
  public:
    modKinaseNoExtrap(double theRate) :
      rate(theRate)
    {}

    double
    getRate(const plx::molInContext& rKinaseContext,
	    const plx::molInContext& rSubstrateContext) const
    {
      return rate;
    }
  };

  // modKinase rate extrapolator that uses reactant masses to extrapolate
  // rates.
  class modKinaseMassExtrap :
    public modKinaseExtrapolator
  {
    double invariant;
  public:
    modKinaseMassExtrap(double rate,
			double leftMolWeight,
			double rightMolWeight) :
      invariant(mzr::bindingInvariant(rate,
				      leftMolWeight,
				      rightMolWeight))
    {}

    double
    getRate(const plx::molInContext& rKinaseContext,
	    const plx::molInContext& rSubstrateContext) const
    {
      plx::cxMol cxKinase(rKinaseContext);
      plx::cxMol cxSubstrate(rSubstrateContext);

      return mzr::bindingRate(invariant,
			      cxKinase.getPlexWeight(),
			      cxSubstrate.getPlexWeight());
    }
  };
}

#endif // MODKINASEEXTRAP_H

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

#ifndef NUCEXEXTRAP_H
#define NUCEXEXTRAP_H

#include "plex/cxMolParam.hh"

namespace nucEx
{
  // This reaction generator actually makes four different kinds of reactions:
  // enabled binding, enabled unbinding, plain binding, and plain unbinding.
  // For this first cut on explicit reaction rate extrapolators, I'm going to
  // assume that the same extrapolator will be used for both kinds of binding
  // reactions and that the same extrapolator will be used for both kinds of
  // unbinding reactions.

  // Base class for (unary) nucleotide unbinding rate extrapolation
  // for nucExRxnGen.
  class nucExUnbindExtrapolator
  {
  public:
    virtual
    ~nucExUnbindExtrapolator(void)
    {}

    virtual double
    getRate(const plx::molInContext& rContext) const = 0;
  };

  // (unary) nucleotide unbinding rate extrapolation for nucExRxnGen.
  class nucExUnbindNoExtrap :
    public nucExUnbindExtrapolator
  {
    double rate;
  public:
    nucExUnbindNoExtrap(double theRate) :
      rate(theRate)
    {}

    double
    getRate(const plx::molInContext& rContext) const
    {
      return rate;
    }
  };

  // Base class for (binary) nucleotide binding rate extrapolation
  // for nucExRxnGen.
  class nucExBindExtrapolator
  {
  public:
    virtual
    ~nucExBindExtrapolator(void)
    {}

    virtual double
      getRate(const plx::molInContext& rContext) const = 0;
  };

  // (binary) nucleotide binding reaction rate extrapolator that
  // does no extrapolation.  Not recommended; included as an example.
  class nucExBindNoExtrap :
    public nucExBindExtrapolator
  {
    double rate;
  public:
    nucExBindNoExtrap(double theRate) :
      rate(theRate)
    {}

    double
    getRate(const plx::molInContext& rContext) const
    {
      return rate;
    }
  };

  // (binary) nucleotide binding reaction rate extrapolator that
  // uses reactant masses to do extrapolation.
  class nucExBindMassExtrap :
    public nucExBindExtrapolator
  {
    double nucleotideMass;
    double invariant;
  public:
    // I notice that, in some other situations such as this, I extracted
    // the masses of the nucleotide and the mol in the constructor.
    nucExBindMassExtrap(double rate,
			double molWeight,
			double nucleotideWeight) :
      nucleotideMass(nucleotideWeight),
      invariant(mzr::bindingInvariant(rate,
				      molWeight,
				      nucleotideMass))
    {}

    double
    getRate(const plx::molInContext& rContext) const
    {
      plx::cxMol cx(rContext);
      
      return mzr::bindingRate(invariant,
			      cx.getPlexWeight(),
			      nucleotideMass);
    }
  };
}

#endif // NUCEXEXTRAP_H

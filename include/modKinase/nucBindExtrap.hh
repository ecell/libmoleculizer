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

#ifndef NUCBINDEXTRAP_H
#define NUCBINDEXTRAP_H

#include "plex/cxMolParam.hh"

namespace kinase
{
  // This header is complicated by the fact that nucBindRxnGen actually
  // generates two different kinds of reactions, both nucleotide binding and
  // unbinding reactions.  These two kinds of reactions need two different
  // rate extrapolators; binding is binary-type but unbinding is unary-type.

  // Base class for (binary) nucleotide binding reaction rate extrapolation.
  class nucBindExtrapolator
  {
  public:
    virtual
    ~nucBindExtrapolator(void)
    {}

    virtual double
    getRate(const plx::molInContext& rContext) const = 0;
  };

  // (binary) nucleotide binding reaction rate extrapolator that does
  // no extrapolation.  Not recommended; included as an example.
  class nucBindNoExtrap :
    public nucBindExtrapolator
  {
    double rate;
  public:
    nucBindNoExtrap(double theRate) :
      rate(theRate)
    {}

    double
    getRate(const plx::molInContext& rContext) const
    {
      return rate;
    }
  };

  // (binary) nucleotide binding reaction rate extrapolator that
  // uses reactant masses to extrapolate binding rate.
  class nucBindMassExtrap :
    public nucBindExtrapolator
  {
    double nucleotideMass;
    double invariant;
  public:
    nucBindMassExtrap(double rate,
		      double nucBindMolWeight,
		      double nucleotideWeight) :
      nucleotideMass(nucleotideWeight),
      invariant(mzr::bindingInvariant(rate,
				      nucBindMolWeight,
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

  // Base class for (unary) nucleotide unbinding reaction rate extrapolation.
  class nucUnbindExtrapolator
  {
  public:
    virtual
    ~nucUnbindExtrapolator(void)
    {}

    virtual double
    getRate(const plx::molInContext& rContext) const = 0;
  };

  // (unary) rate extrapolation for nucleotide unbinding.
  class nucUnbindNoExtrap :
    public nucUnbindExtrapolator
  {
    double rate;
  public:
    nucUnbindNoExtrap(double theRate) :
      rate(theRate)
    {}

    double
    getRate(const plx::molInContext& rContext) const
    {
      return rate;
    }
  };
}

#endif // NUCBINDEXTRAP_H

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

#ifndef CFT_UNIMOLEXTRAP_H
#define CFT_UNIMOLEXTRAP_H

#include "cml/cptModMol.hh"
#include "clx/cptPlexSpecies.hh"
#include "clx/cptPlexFamily.hh"

namespace cft
{
  class uniMolExtrapolator
  {
  public:
    virtual
    ~uniMolExtrapolator(void)
    {}

    virtual
    double
    getRate(const cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily>& rContext) const = 0;
  };

  class uniMolNoExtrap :
    public uniMolExtrapolator
  {
    double rate;
  public:
    uniMolNoExtrap(double theRate) :
      rate(theRate)
    {}

    double
    getRate(const cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily>& rContext) const
    {
      return rate;
    }
  };

  class uniMolMassExtrap :
    public uniMolExtrapolator
  {
    // Pointer to "massive" part of auxiliary reactant species; null if there
    // is no auxiliary reactant species.
    const fnd::massive* pMassive;

    // The rate if the reaction is unary (0 == pMassive) but the binding
    // invariant if the reaction is binary.
    double rateOrInvariant;

  public:
    // For creating unary reactions.
    uniMolMassExtrap(double theRate) :
      pMassive(0),
      rateOrInvariant(theRate)
    {}

    // For creating binary reactions.
    uniMolMassExtrap(double theRate,
		     const cml::cptModMol* pEnablingMol,
		     const fnd::massive* pMassiveAuxiliarySpecies);

    double
    getRate(const cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily>& rContext) const;
  };
}

#endif //  CFT_UNIMOLEXTRAP_H

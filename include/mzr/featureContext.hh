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

#ifndef FEATURECONTEXT_H
#define FEATURECONTEXT_H

/*! \file featureContext.hh
  \ingroup featureGroup
  \brief Defines featureContext, an entry in a feature's contexts vector. */

// For the definition of pair.
#include <utility>

namespace mzr
{
  /*! \ingroup featureGroup
    \brief An entry in the contexts vector of a feature.

    A member of this class points to a species and tells how a feature
    is presented by that species.  For example, it might give a species
    of complex and say that we're looking at the third binding in that
    complex.

    Features keep a vector of all of their contexts, so that
    reactionFamilies that are sensitive to the feature can generate new
    reactions for each new presentation of the feature.  */
  template<class bearerSpecies, class featureSpec>
  class featureContext :
    public std::pair<bearerSpecies*, featureSpec>
  {
  public:


    /*!  \brief Abstract base class for reaction generator.
    
    This is the means by which a reactionFamily attaches to a feature.
    Each rxnGen of a reactionFamily will belong to a unique class
    derived from this.  That is, each rxnGen of each reactionFamily
    will have to give an implementation of makeReactions. */
    class rxnGen
    {
    public:
      virtual ~rxnGen(void)
      {}

      virtual void
      makeReactions(const featureContext& rContext) const = 0;
    };

    featureContext(void)
    {}
    
    featureContext(bearerSpecies* pSpecies,
		   const featureSpec& rSpec) :
      std::pair<bearerSpecies*, featureSpec>(pSpecies, rSpec)
    {}

    //! Returns a pointer to the context species.
    bearerSpecies*
    getSpecies(void) const
    {
      return this->first;
    }

    /*! \brief Returns the way that the species displays the feature.

    For example, if the feature is a binding, then this returns the
    index of the binding in the plexSpecies. */
    const featureSpec&
    getSpec(void) const
    {
      return this->second;
    }
  };
}

#endif

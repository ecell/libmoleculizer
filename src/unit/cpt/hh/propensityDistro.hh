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

#ifndef PROPENSITYDISTRO_H
#define PROPENSITYDISTRO_H

#include <list>
#include "utl/xcpt.hh"
#include "utl/gsl.hh"

namespace cpt
{
  class cptReaction;
  
  class propensityDistro :
      public std::list<std::pair<double, cptReaction*> >
  {
    utl::gsl::uniformSampler uniform;
    
    long double totalPropensity;

    // This is for debugging with sumPropensity below.
    class addPropensity :
      public std::unary_function<std::pair<double, cptReaction*>, void>
    {
      long double& rSum;
    public:
      addPropensity(long double& rSumOfPropensities) :
	rSum(rSumOfPropensities)
      {
      }

      void
      operator()(const argument_type& rPropensityReactionPair) const
      {
	rSum += rPropensityReactionPair.first;
      }
    };

  public:
    propensityDistro(utl::gsl::autoGslRng& rAutoGslRng) :
      uniform(rAutoGslRng),
      totalPropensity(0.0)
    {}

    long double
    getTotalPropensity(void)
    {
      return totalPropensity;
    }

    // For debugging local modification of total propensity.
    // Check this sum against total propensity.
    long double
    sumPropensity(void)
    {
      long double theSum = 0.0;
      std::for_each(begin(),
		    end(),
		    addPropensity(theSum));
      return theSum;
    }

    // This is so that a reaction, by looking at its former propensity
    // and its new propensity, can update the total propensity.
    void
    updateTotalPropensity(double delta)
    {
      totalPropensity += delta;
    }
    
    // Move the reaction/propensity with the given iterator in the
    // propensity distribution to the front of the propensity distribution.
    //
    // This is supposed to be a cheap way of keeping fast reactions toward
    // the front of the propensity distribution without actually doing any
    // kind of sorting.
    iterator
    moveToFront(iterator iDistroEntry)
    {
      push_front(*iDistroEntry);
      erase(iDistroEntry);
      return begin();
    }

    iterator
    sample(void)
      throw(utl::xcpt);

    // Adds the reaction with the assumption that it has
    // propensity 0.  This should be done before the reaction
    // has ever "responded," so that the propensity of 0 with which
    // it is entered here agrees with the reaction's lastCalculatedPropensity.
    void
    addReaction(cptReaction* pCptReaction);
  };
}
  
#endif // PROPENSITYDISTRO_H

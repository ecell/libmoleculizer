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

#ifndef FEATUREMAP_H
#define FEATUREMAP_H

/*! \file featureMap.hh
  \ingroup featureGroup
  \brief Defines container for features with common spec type. */

#include <map>
#include "mzr/feature.hh"

namespace mzr
{
/*! \ingroup featureGroup
  \brief A container for features with common spec type.

  A species family uses this class to group together all its
  features for sites, or all its features for bindings, etc.

  This class is a leftover from a previous implementation in which
  features had a more dynamic role and had more member functions to be
  called by a connected speciesFamily.  It does no harm, but I
  probably wouldn't have defined it just to encapsulate such limited
  functionality.

  \todo Might want to eliminate this nearly un-called-for class.
*/
template<class bearerSpecies, class featureSpec>
class featureMap :
  public std::map<featureSpec, feature<bearerSpecies, featureSpec>*>
{
public:

  typedef
  typename std::map<featureSpec, feature<bearerSpecies, featureSpec>*>::value_type
  value_type;

private:
  
  //! Function class in support of notifyNew.
  class notifyFeature :
    std::unary_function<featureMap::value_type, void>
  {
    bearerSpecies* pSpecies;
  public:
    notifyFeature(bearerSpecies* pNewSpecies) :
      pSpecies(pNewSpecies)
    {}
    
    void operator()(const typename featureMap::value_type& rEntry)
    {
      feature<bearerSpecies, featureSpec>* pFeature = rEntry.second;
      const featureSpec& rSpec = rEntry.first;

      pFeature->notifyNew(pSpecies, rSpec);
    }
  };
  
public:

  //! Add a feature to be notified when a new species appears.
  void
  addFeature(const featureSpec& rSpec,
	     feature<bearerSpecies, featureSpec>* pFeature)
  {
    bool insertSucceeded
      = insert(featureMap::value_type(rSpec, pFeature)).second;
    
    if(! insertSucceeded)
      throw featureAlreadyMappedXcpt();
  }

  //! Notify each feature targeted by the map of the new species.
  void
  notifyNew(bearerSpecies* pNewSpecies) const
  {
    for_each(begin(),
	     end(),
	     notifyFeature(pNewSpecies));
  }
};
}

#endif

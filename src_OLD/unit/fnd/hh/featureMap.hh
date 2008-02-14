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
#include "utl/xcpt.hh"
#include "fnd/feature.hh"

namespace fnd
{
  // It would probably be better to template on the feature type,
  // since the feature is itself templated on the context type.
  template<class contextT>
  class featureMap :
    public std::map<typename contextT::contextSpec, feature<contextT>* >,
    public sensitive<newSpeciesStimulus<typename contextT::speciesType> >
  {
    class notifyFeature :
      std::unary_function<typename featureMap::value_type, void>
    {
      const typename featureMap::stimulusType& rStimulus;
      
    public:
      notifyFeature(const typename featureMap::stimulusType& rNewSpeciesStimulus) :
	rStimulus(rNewSpeciesStimulus)
      {}
    
      void operator()(const typename featureMap::value_type& rEntry) const
      {
	const typename contextT::contextSpec& rSpec = rEntry.first;
	feature<contextT>* pFeature = rEntry.second;

	contextT newContext(rStimulus.getSpecies(),
			    rSpec);

	newContextStimulus<contextT> stim(newContext,
					  rStimulus.getNotificationDepth());
	
	pFeature->respond(stim);
      }
    };

  public:
    typedef newSpeciesStimulus<typename contextT::speciesType> stimulusType;
    
    //! Add a feature to be notified when a new species appears.
    bool
    addFeature(const typename contextT::contextSpec& rSpec,
	       feature<contextT>* pFeature)
      throw(utl::xcpt)
    {
      return insert(typename featureMap::value_type(rSpec, pFeature)).second;
    }

    //! Notify each feature targeted by the map of the new species.
    void
    respond(const typename featureMap::stimulusType& rStimulus)
    {
      for_each(this->begin(),
	       this->end(),
	       notifyFeature(rStimulus));
    }
  };
}

#endif

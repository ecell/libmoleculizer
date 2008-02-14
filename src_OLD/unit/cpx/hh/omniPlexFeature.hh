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

#ifndef OMNIPLEXFEATURE_H
#define OMNIPLEXFEATURE_H

/*! \file omniPlexFeature.hh
  \ingroup omniGroup
  \brief Defines subcomplex feature. */

#include "fnd/feature.hh"
#include "fnd/multiSpeciesDumpable.hh"
#include "fnd/newContextStimulus.hh"
#include "cpx/cxOmni.hh"

namespace cpx
{
  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  class omniPlexFeature :
    public fnd::feature<cxOmni<molT,
			       plexSpeciesT,
			       plexFamilyT,
			       omniPlexT> >
  {
  public:
    typedef plexSpeciesT plexSpeciesType;
    typedef plexFamilyT plexFamilyType;
    typedef omniPlexT omniPlexType;

    typedef typename plexSpeciesT::msDumpableType dumpableType;

    typedef andPlexQueries<plexSpeciesType,
			   omniPlexType> stateQueryType;

    typedef
    typename cpx::cxOmni<molT, plexSpeciesT, plexFamilyT, omniPlexT>
    contextType;

    typedef fnd::newContextStimulus<contextType> stimulusType;

  private:
    // If a dumpable is really attached to this feature, then
    // this points to it; otherwise null.
    dumpableType* pDumpable;
    
  public:
    omniPlexFeature(void) :
      pDumpable(0)
    {}

    // Generate reactions and add new species to the dumpable,
    // if any.
    //
    // Overrides fnd::feature<cpx::cxOmni>::respond to add new species
    // to possible dumpable.
    virtual
    void
    respond(const typename omniPlexFeature::stimulusType& rNewFeatureContext);

    // To "turn on" dumping of the species in this omniplex.
    // These aren't query-based dumpables, since the omniplex
    // itself does all the querying.
    void
    setDumpable(typename omniPlexFeature::dumpableType* ptrDumpable)
    {
      pDumpable = ptrDumpable;
    }
  };

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  void
  omniPlexFeature<molT,
		  plexSpeciesT,
		  plexFamilyT,
		  omniPlexT>::
  respond(const typename omniPlexFeature::stimulusType& rStim)
  {
    const typename omniPlexFeature::contextType& rNewContext
      = rStim.getContext();
    
    // Does the new species satisfy the omni's state query?
    omniPlexType* pOmni = rNewContext.getOmni();
    const typename omniPlexFeature::stateQueryType& rQuery
      = *(pOmni->getStateQuery());

    plexSpeciesType* pSpecies
      = rNewContext.getSpecies();
    const subPlexSpec<omniPlexType>& rSpec
      = rNewContext.getSpec();
    if(rQuery.applyTracked(*pSpecies,
			   rSpec))
      {
	// Notify reaction generators 
	fnd::feature<contextType>::respond(rStim);

	// Notify dumpable if any.
	if(pDumpable)
	  {
	    // Note that this just gets back the newSpeciesStimulus.
	    // This really stinks.
	    fnd::newSpeciesStimulus<plexSpeciesType>
	      dumpStim(pSpecies,
		       rStim.getNotificationDepth());
	    
	    pDumpable->respond(dumpStim);
	  }
      }
  }
}

#endif

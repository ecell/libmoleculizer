//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

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

    typedef andPlexQueries<plexSpeciesType,
    omniPlexType> stateQueryType;

    typedef
    typename cpx::cxOmni<molT, plexSpeciesT, plexFamilyT, omniPlexT>
    contextType;

    typedef fnd::newContextStimulus<contextType> stimulusType;

private:

public:
    omniPlexFeature( void )
    {}

    // Generate reactions
    //
    // Overrides fnd::feature<cpx::cxOmni>::respond to add new species
    // to possible dumpable.
    virtual
    void
    respond( const typename omniPlexFeature::stimulusType& rNewFeatureContext );
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
respond( const typename omniPlexFeature::stimulusType& rStim )
{
    const typename omniPlexFeature::contextType& rNewContext
    = rStim.getContext();

    // Does the new species satisfy the omni's state query?
    omniPlexType* pOmni = rNewContext.getOmni();
    const typename omniPlexFeature::stateQueryType& rQuery
    = * ( pOmni->getStateQuery() );

    plexSpeciesType* pSpecies
    = rNewContext.getSpecies();
    const subPlexSpec<omniPlexType>& rSpec
    = rNewContext.getSpec();
    if ( rQuery.applyTracked( *pSpecies,
                              rSpec ) )
    {
        // Notify reaction generators
        fnd::feature<contextType>::respond( rStim );
    }
}
}

#endif

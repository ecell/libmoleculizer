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

#ifndef PLEXFEATURE_H
#define PLEXFEATURE_H

/*! \defgroup plexFeatureGroup Features
  \ingroup plexGroup

  \brief Features notifying reaction families when new complex species
  appear. */

/*! \file plexFeature.hh
  \ingroup plexFeatureGroup
  \brief Typedefs of features and "species in context" for complexes. */

#include "mzr/feature.hh"
#include "plex/plexSpec.hh"
#include "plex/subPlexSpec.hh"
#include "plex/plexSpecies.hh"

namespace plx
{
  /*! \name Feature typedefs
    \ingroup plexFeatureGroup
    \brief Features of complexes.

    These objects notify reaction families when a new species of complex
    that "has the feature" appears.  The reaction family is notified
    with a "feature in context," the context being a species and a
    "spec" that tells how the species "presents" the feature.

    For example, suppose a new species S of complex appears that has
    binding between the third binding site on X and the first binding
    site on Y.  The structural family to which the new complex belongs
    notifies the feature for (X, site 3) <--> (Y, site 1) bindings of
    the new species.  (Note that many other structural families of
    complexes will also notify this feature from time to time.)  The
    feature then notifies the decomposition reaction family for this
    binding that an S has appeared, and that the binding in question is
    the 10th binding in S.  This lets the decomposition reaction family
    construct the reaction decomposing S's at the 10th binding, using
    what the reaction family already knows about binding of the 3rd site
    on X to the 1st site on Y. Conceivably, other reactions might be
    interested in this binding and would also be attached to the (X,
    site 3) <--> (Y, site 1) binding feature to hear when new species
    appear that contain this particular binding. */
  //@{
  typedef mzr::feature<plexSpecies, plexSiteSpec> siteFeature;
  typedef mzr::feature<plexSpecies, plexBindingSpec> bindingFeature;
  typedef mzr::feature<plexSpecies, plexMolSpec> molFeature;
  typedef mzr::feature<plexSpecies, subPlexSpec> omniFeature;
  //@}  

  /*! \name Feature context typedefs
    \ingroup plexFeatureGroup
    \brief Features in context for complexes.

    These simple objects tell a reaction family how a particular feature
    is "presented" by a species of complex.  For example, a "binding
    in context" is a species of complex, together with the index
    of the binding (in the plexFamily's paradigm complex.). */
  //@{
  typedef mzr::featureContext<plexSpecies, plexSiteSpec> siteInContext;
  typedef mzr::featureContext<plexSpecies, plexBindingSpec> bindingInContext;
  typedef mzr::featureContext<plexSpecies, plexMolSpec> molInContext;
  typedef mzr::featureContext<plexSpecies, subPlexSpec> omniInContext;
  //@}
}

#endif

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

#ifndef FEATURE_H
#define FEATURE_H

/*! \defgroup featureGroup Features
  \ingroup mzrGroup
  \brief Template classes to notify reaction families when new species appear.

  A feature is (associated to) a characteristic of a molecule that
  might make that molecule a substrate of a reaction.  For example, a
  free binding site S is associated to a feature.  Each structural
  family of complexes with a free S site connects to S's feature.
  Each family of reactions that binds S to something else also
  connects to S's feature.  Then, when a new complex with a free S
  site appears, the reaction family is notified and the reactions that
  bind the new complex to other complexes at S can be created.

  The feature is the means by which any one of many species families
  can notify all of a number of reaction families when the species family
  has a new member species.
  
  In a complex, the mols, bindings, free binding sites, and any
  recognized subcomplexes are associated with features. */

/*! \file feature.hh
  \ingroup featureGroup
  \brief Defines template class feature. */

#include <vector>
#include <functional>
#include <algorithm>
#include "mzr/featureContext.hh"
#include "mzr/util.hh"

namespace mzr
{
  /*! \ingroup featureGroup
    \brief Notifies reactions of new species.

    A feature is some aspect of a species that is the subject of a
    reactionFamily.  It is the reason that the reactionFamily is
    interested in a speciesFamily.

    For example, the free binding sites on a complex are features.  They
    are the reason that the family of dimerization reactions is
    interested in a structural family of complexes, which all share the
    same free binding sites. Many plexFamilies will connect to the same
    free binding site feature, since many structural families will have
    the same binding site free.  Several reaction families might connect
    to the free binding site feature if the binding site can bind with
    several partners.

    Subcomplexes are features, since some families of reactions are
    interested in all complexes that contain some particular subcomplex. */
  template<class bearerSpecies, class featureSpec>
  class feature
  {
  public:

    typedef featureContext<bearerSpecies, featureSpec> context;
    typedef typename context::rxnGen rxnGen;
    typedef typename std::vector<context>::iterator contextIterator;

  private:
  
    class notifyGenerator :
      std::unary_function<rxnGen*, void>
    {
      const context& rContext;
    public:
      notifyGenerator(const context& rNewContext) :
	rContext(rNewContext)
      {}

      void operator()(rxnGen* pRxnGen) const
      {
	pRxnGen->makeReactions(rContext);
      }
    };

  public:
    virtual ~feature(void)
    {}

    /*! \brief Species that have the feature and how they express it.

    For example, if the feature is a binding, then each context gives
    the plexSpecies and the index of the binding in the
    plexSpecies.

    This member is public because reaction generators (rxnGen) have to
    be able to run through all the species that are connected to a
    feature.  For example, dimerization needs to know all the possible
    left-hand binding partner species for a new species displaying the
    dimerization's right-hand free site feature.  It finds these
    possible binding partners in the left-hand feature's contexts
    vector.  */
    std::vector<context> contexts;

    /*! \brief Parties interested in new species bearing the feature.

    Each reactionFamily that might need to generate reactions has
    reaction generators, rxnGen, that will generate the needed reactions
    when notified by the feature.  When the reactionFamily is constructed,
    it is "hooked up" to the features that it's interested in by pushing
    its rxnGen reaction generators onto this vector in the features.

    When notified by one of its speciesFamilies of a new species, the
    feature runs down its vector of rxnGens, notifying them in turn
    about the new species.

    Note that reaction generators cannot be memory managed in the
    obvious way by the features that they are connected to, since a
    reaction generator can be connected to more than one feature.
    Moreover, the reaction generators associated with binary reactions
    are not independent objects; they are consituents of a binary
    reaction generator pair.  */
    std::vector<rxnGen*> rxnGens;

    /*! \brief A speciesFamily uses this to notify feature of new species.

    This is virtual because omniplex features need to notify their
    associated dumpables of new species. */
    virtual void
    notifyNew(bearerSpecies* pNewSpecies,
	      const featureSpec& rSpec)
    {
      context newContext(pNewSpecies,
			 rSpec);
    
      // Add the new context to the vector of all contexts in which
      // this feature occurs.
      contexts.push_back(newContext);

      // Have every rxnGen generate all the new reactions in its
      // associated reaction family that involve the new species. Each
      // new reaction is added to the sensitivity list of the new species,
      // and to the sensitivity lists of all the other species that
      // are involved.
      for_each(rxnGens.begin(),
	       rxnGens.end(),
	       notifyGenerator(newContext));
    }

    /*! \brief A reactionFamily uses this to connect to feature.

    A reaction family may have more than one reaction generator.  For
    example, the dimerization family, dimerizeFam, has two (bundled
    together for convenience.)  One is connected to the left-hand
    free site (feature) of the dimerization (s) and the other is connected
    to the right-hand free site (feature) with this function. */
    void
    addRxnGen(rxnGen* pRxnGen)
    {
      rxnGens.push_back(pRxnGen);
    }
  };
}

#endif

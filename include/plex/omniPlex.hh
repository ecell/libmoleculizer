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

#ifndef OMNIPLEX_H
#define OMNIPLEX_H


#include "plex/omniPlexFeature.hh"
#include "plex/alloSiteQuery.hh"
#include "plex/omniStructureQuery.hh"

namespace plx
{
  class omniPlex
  {
    // Pointer to the parent node under which this omniPlex
    // was parsed.  This is here to make it possible for
    // modules to get to the C++ omniplexes that are parsed
    // from nodes that the modules point out to the plex module.
    //
    // Formerly, one only needed the plexFamily associated
    // to a particular omniplex construct in the input file,
    // and this was obtained by "unification," i.e. by looking the
    // complex up in the structural database.  Now, in addition
    // to doing that, module code will need to run down all the omniPlexes
    // associated to the plexFamily, looking for the one that was
    // parsed for module by the plexUnit.
    xmlpp::Node* pNode;

    // Pointer to the plexFamily of this omniplex.
    //
    // One side-effect of this way of working is that one can't get from a
    // plexFamily to all of its associated omniplex stuff, since the same
    // plexFamily may be an omniplex in more than one way now.
    plexFamily* pFamily;

    // Reactions that want to listen for species that satisfy the
    // structural queries should listen to this feature.
    omniPlexFeature subPlexFeature;

    // Structural queries that must be statisfied for this
    // omniplex to be connected to a plexFamily as one of its features.
    //
    // In order for compound queries like this to be copyable, it's necessary
    // (really easiest) to memory-manage queries centrally, in the plexUnit.
    // This andOmniStructureQueries object should be registered with plexUnit
    // before incorporation itno the omniplex.
    andOmniStructureQueries* pStructureQueries;

    // State queries that must be satisfied for a plexSpecies
    // to be (allo)sterically modified as given by the siteShapeMap.
    //
    // In order for compound queries like this to be copyable, it's necessary
    // (really easiest) to memory-manage queries centrally, in the plexUnit.
    // This andPlexQueries object should be registered with plexUnit before
    // incorporation into the omniPlex.
    andPlexQueries* pStateQueries;

    // Map giving allosteric site shapes.  When a plexFamily is recognized, it
    // may pass the structureQueries above and be connected to this omni.
    // Then, when a new plexSpecies in that family is created, its state may
    // pass the stateQueries above.  If so, then the shapes of the sites in
    // the new plexSpecies that correspond to the sites in this map are set to
    // the shapes given by this map.
    //
    // This is only used in allosteric omnis; it remains empty as constructed
    // otherwise.
    siteToShapeMap alloSiteMap;
    
  public:
    omniPlex(xmlpp::Node* pParentNode,
	     plexFamily* pPlexFamily,
	     andOmniStructureQueries* pOmniStructureQueries,
	     andPlexQueries* pOmniStateQueries);

    // This node pointer becomes useless after the parsing phase,
    // when the DOM object is deleted.
    xmlpp::Node*
    getParentNode(void) const
    {
      return pNode;
    }

    plexFamily*
    getFamily(void) const
    {
      return pFamily;
    }

    // For use in plexFamily::connectToFeatures.
    const andOmniStructureQueries&
    getStructureQuery(void) const
    {
      return *pStructureQueries;
    }

    // Used in omniPlexFeature::notifyNew.  Perhaps omniPlexFeature should
    // be a friend?
    const andPlexQueries*
    getStateQuery(void) const
    {
      return pStateQueries;
    }

    // Replaces plexFamily::getSubPlexFeature, now that a complex can be omni
    // in more than one way.
    omniPlexFeature*
    getSubPlexFeature(void)
    {
      return & subPlexFeature;
    }

    // For use in plex allostery routines.
    const siteToShapeMap&
    getSiteToShapeMap(void) const
    {
      return alloSiteMap;
    }

    // For use in the parser.
    siteToShapeMap&
    getSiteToShapeMap(void)
    {
      return alloSiteMap;
    }
  };
}

#endif // OMNIPLEX_H

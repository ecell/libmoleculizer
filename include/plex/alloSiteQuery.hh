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

#ifndef ALLOSITEQUERY_H
#define ALLOSITEQUERY_H

// #include "plex/plexQuery.hh"

#include <map>
#include "mol/siteShape.hh"
#include "plex/plexXcpt.hh"

// This file defines classes used in plexFamily to record the allosteric
// shapes of binding sites when the mols in the plex satisfy different
// conditions on their states.  This is to help enable allosteric-omnis to
// include state specifications.

namespace plx
{
  class plexSiteSpec;
  class subPlexSpec;
  class plexParam;
  class plexIsoPair;
  class andPlexQueries;

  // This class is used to record the allosteric shapes of the binding
  // sites in a complex.
  class siteToShapeMap :
    public std::map<plexSiteSpec, bnd::siteParam>
  {
  public:
    bnd::siteParam
    getSiteShape(const plexSiteSpec& rPlexSiteSpec) const
      throw(unmappedSiteSpecXcpt);

    // Forces insertion of the one given site-shape pair.
    void
    setSiteShape(const plexSiteSpec& rPlexSiteSpec,
		 const bnd::siteParam& rSiteParam);

    // Same as above, but does the insertion through the
    // given subPlexSpec.
    void
    setSiteShape(const plexSiteSpec& rPlexSiteSpec,
		 const bnd::siteParam& rSiteParam,
		 const subPlexSpec& rSubPlexSpec);

    // Decided, for reasons unknown, to use overloading for the
    // version of this function that works through the injection
    // of a subcomplex into a complex.

    // Forces insertion of all the site-to-shape entries from the
    // given rSiteToShapeMap.
    void
    setSiteShapes(const siteToShapeMap& rSiteToShapeMap);

    // Forces insertion of all the site-to-shape entries from the given
    // rSiteToShapeMap, but pushes all of the site specifications
    // through the given injection.
    void
    setSiteShapes(const siteToShapeMap& rSiteToShapeMap,
		  const subPlexSpec& rSubPlexSpec);
  };

  // This lists the allosteric shapes that a complex takes on when
  // its components satisfy state queries.  This optimizes traversal,
  // since we must run through all the queries for every new species
  // in the plex family.
  class queryAllosteryList :
    public std::list<std::pair<const andPlexQueries*, siteToShapeMap> >
  {
  public:
    // This is used to add allostery to a plexFamily.
    // 
    // One adds a plexQuery and a specification of all the shapes of the
    // complex's binding sites.  When a new species of complex in the
    // plexFamily appears, then then the query is applied.  If it tests true,
    // then the associated shapes are overlaid on the new species's binding
    // site shapes.
    void
    addQueryAndMap(const andPlexQueries* pQuery,
		   const siteToShapeMap& rSiteToShapeMap);
    
    // Decided, for reasons unknown, to use overloading for the
    // version of this function that works through the injection
    // of a subcomplex into a complex.

    // Runs down the list, seeing which queries are satisfied
    // by the plexParam.  When a query, for now a query about
    // the states of the mols in the plex, is satisfied, then
    // the site shapes in the plexParam are set from the siteShapeMap
    // associated to the query.
    //
    // This is used in the construction of new plexSpecies that belong to a
    // plexFamily for which allosteric states have been declared (an
    // allosteric-plex.)
    void
    setSatisfiedQuerySiteShapes(plexParam& rPlexParam) const;

    // Same as the above, but the queries are pushed forward through
    // the injection, as are the changes in site shapes.
    //
    // This is used in the construction of new plexSpecies in which
    // allosteric subcomplexes have been detected.
    void
    setSatisfiedQuerySiteShapes(plexParam& rParam,
				const subPlexSpec& rSubPlexSpec) const;
  };
}

#endif // ALLOSITEQUERY_H

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

#include "plex/omniStructureQuery.hh"

namespace plx
{
  class applyNotStructureQuery :
    public std::unary_function<omniStructureQuery*, bool>
    {
      const plex& rPlex;
      const plexIsoPair& rInj;
    public:
      applyNotStructureQuery(const plex& rArgumentPlex,
			     const plexIsoPair& rInjection) :
	rPlex(rArgumentPlex),
	rInj(rInjection)
      {}

      bool
      operator()(const omniStructureQuery* pQuery) const
      {
	const omniStructureQuery& rQuery = *pQuery;
	return ! rQuery(rPlex,
			rInj);
      }
    };

  // AND's the values of structural queries.
  bool
  andOmniStructureQueries::
  operator()(const plex& rPlex,
	     const plexIsoPair& rInjection) const
  {
    return (queries.end() == find_if(queries.begin(),
				     queries.end(),
				     applyNotStructureQuery(rPlex,
							    rInjection)));
  }

  class bindingContainsSite :
    public std::unary_function<plexBinding, bool>
  {
    // The site being searched for among the bindings
    // AFTER it has been translated by the omni imbedding.
    const plexSiteSpec& rQuerySiteSpec;

  public:
    bindingContainsSite(const plexSiteSpec& rSiteSpec) :
      rQuerySiteSpec(rSiteSpec)
    {}
    
    bool
    operator()(const plexBinding& rBinding) const
    {
      return (rBinding.leftSite() == rQuerySiteSpec)
	|| (rBinding.rightSite() == rQuerySiteSpec);
    }
  };

  bool
  omniFreeSiteQuery::
  operator()(const plex& rPlex,
	     const plexIsoPair& rInjection) const
  {
    // Translate the spec of the site we want to check into the indexing of
    // the given plex using the injection of the omniplex into the given plex
    // where it has been recognized.
    plexSiteSpec translatedSiteSpec
      = rInjection.forward.applyToSiteSpec(freeSiteSpec);

    // Return true if no binding in the plex contains the translated site
    // spec; i.e. if the site specified by freeSiteSpec in the omni
    // is also free in the plex where the omni is embedded.
    return
      (rPlex.bindings.end()
       == std::find_if(rPlex.bindings.begin(),
		       rPlex.bindings.end(),
		       bindingContainsSite(translatedSiteSpec)));
  }
}


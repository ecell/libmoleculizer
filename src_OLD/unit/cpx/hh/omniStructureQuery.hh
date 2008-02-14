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

#ifndef OMNISTRUCTUREQUERY_H
#define OMNISTRUCTUREQUERY_H

#include "cpx/plexIso.hh"
#include "cpx/ftrSpec.hh"

namespace cpx
{
  template<class plexT>
  class omniStructureQueryArg
  {
  public:
    const plexT& rPlex;
    const plexIso& rInjection;

    omniStructureQueryArg(const plexT& rStructurePlex,
			  const plexIso& rInjectionIso) :
      rPlex(rStructurePlex),
      rInjection(rInjectionIso)
    {}
  };

  // Examines matching plex to see if specified free sites in the
  // omni are also free sites in the matching plex.
  template<class plexT>
  class omniFreeSiteQuery :
    public fnd::query<omniStructureQueryArg<plexT> >
  {
    siteSpec freeSiteSpec;

  public:
    omniFreeSiteQuery(const siteSpec& rFreeSiteSpec) :
      freeSiteSpec(rFreeSiteSpec)
    {}

    ~omniFreeSiteQuery(void)
    {}
    
    bool
    operator()(const omniStructureQueryArg<plexT>& rArg) const;
  };
}

#include "cpx/omniStructureQueryImpl.hh"

#endif // OMNISTRUCTUREQUERY_H

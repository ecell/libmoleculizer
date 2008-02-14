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

#ifndef CPX_HASHMOLREC_H
#define CPX_HASHMOLREC_H

#include "utl/linearHash.hh"

namespace cpx
{
  // Function class to hash plexes by "topological sorting."
  template<class molT>
  class hashMolRec
  {
    const basicPlex<molT>& rPlex;
    const typename std::map<siteSpec, int>& rSiteToBindings;
    typename std::set<int>& rMolsSeen;
  
  public:
    hashMolRec(const basicPlex<molT>& refPlex,
	       const typename std::map<siteSpec, int>& refSiteToBindings,
	       typename std::set<int>& refMolsSeen) :
      rPlex(refPlex),
      rSiteToBindings(refSiteToBindings),
      rMolsSeen(refMolsSeen)
    {}
  
    size_t operator()(int molNdx,
		      int depth) const;
  };
}

#include "cpx/hashMolRecImpl.hh"

#endif // CPX_HASHMOLREC_H

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

#ifndef MOL_UNKSITESHAPEXCPT_H
#define MOL_UNKSITESHAPEXCPT_H

#include "utl/xcpt.hh"
#include "utl/dom.hh"
#include "mol/mzrBndSite.hh"
#include "mol/mzrMol.hh"

namespace bnd
{
  // Exception thrown when the user refers to a binding site shape
  // by an incorrect name.
  class unkSiteShapeXcpt : 
    public utl::xcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const mzrBndSite& rBindingSite,
	  const mzrMol* pMol,
	  const std::string& rBadSiteShapeName);
  public:
    unkSiteShapeXcpt(const xmlpp::Node* pOffendingNode,
			 const mzrBndSite& rBindingSite,
			 const mzrMol* pMol,
			 const std::string& rBadSiteShapeName) :
      utl::xcpt(mkMsg(pOffendingNode,
		      rBindingSite,
		      pMol,
		      rBadSiteShapeName))
    {}
  };
}

#endif // MOL_UNKSITESHAPEXCPT_H

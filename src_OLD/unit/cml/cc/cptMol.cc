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

#include <sstream>
#include "cml/cptMol.hh"
#include "cml/unkSiteXcpt.hh"
#include "clx/cptPlexFamily.hh"

namespace cml
{
  cptMol::
  cptMol(const std::string& rName,
	 const std::vector<cptBndSite>& rSites,
	 const std::vector<double>& rBoundaryRates) :
    cpx::basicMol<cptBndSite>(rName,
			      rSites),
    boundaryRates(rBoundaryRates)
  {}

  int
  cptMol::
  mustFindSite(const std::string& rSiteName,
	       xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    int siteNdx = -1;

    if(! findSite(rSiteName,
		  siteNdx))
      throw(unkSiteXcpt(rSiteName,
			pRequestingNode));

    return siteNdx;
  }

  cptBndSite*
  cptMol::
  mustGetSite(const std::string& rSiteName,
	      xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    cptBndSite* pSite = getSite(rSiteName);

    if(! pSite)
      throw(unkSiteXcpt(rSiteName,
			pRequestingNode));

    return pSite;
  }

  std::string
  cptMol::
  genInstanceName(int molInstanceNdx) const
  {
    std::ostringstream oss;
    oss << "mzr-mol_"
	<< molInstanceNdx;
    return oss.str();
  }
}

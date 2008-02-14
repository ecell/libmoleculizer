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

#ifndef CPX_NOKINETICCONSTSXCPT_H
#define CPX_NOKINETICCONSTSXCPT_H

#include <utl/xcpt.hh>

namespace cpx
{
  // Thrown in plexAllostery.cc,
  // plexFamily::insertBindingParam::operator(),when the binding param of each
  // binding in a new plex species is looked up.
  //
  // The first form of the exception gives a not-very-informative message,
  // but it isn't very likely (impossible?) now that kinetics for pairs
  // of binding sites mentioned in dimerization-gen's are defaulted for
  // all possible pairs of site shapes.
  //
  // Thrown in second form in plexConnect.cc, plexFamily::behaviorizeBinding,
  // when a new plexFamily is connected with each of its binding features.
  // This is likely to happen regularly, when the user fails to give a
  // dimerization-gen for some binding in a complex that she specifies.
  //
  // The common thread is that the binding features and kinetics are
  // both stored in databases in the dimerUnit; these databases are
  // filled from dimerization-gen constructs.
  class noKineticConstsXcpt : 
    public utl::xcpt
  {
    static std::string
    mkShapesOnlyMsg(const std::string& rLeftSiteShapeName,
		    const std::string& rRightSiteShapeName);

    static std::string
    mkMolsAndSitesMsg(const std::string& rLeftMolName,
		      const std::string& rLeftSiteName,
		      const std::string& rRightMolName,
		      const std::string& rRightSiteName);

    static std::string
    mkFullMsg(const std::string& rLeftMolName,
	      const std::string& rLeftSiteName,
	      const std::string& rLeftSiteShapeName,
	      const std::string& rRightMolName,
	      const std::string& rRightSiteName,
	      const std::string& rRightSiteShapeName);

    noKineticConstsXcpt(const std::string& rMsg) :
      utl::xcpt(rMsg)
    {}
      
  public:
    static
    noKineticConstsXcpt
    shapesOnly(const std::string& rLeftSiteShapeName,
	       const std::string& rRightSiteShapeName)
      {
	return
	  noKineticConstsXcpt
	  (mkShapesOnlyMsg(rLeftSiteShapeName,
			 rRightSiteShapeName));
      }

    static
    noKineticConstsXcpt
    molsAndSites(const std::string& rLeftMolName,
		 const std::string& rLeftSiteName,
		 const std::string& rRightMolName,
		 const std::string& rRightSiteName)
    {
      return
	noKineticConstsXcpt
	(mkMolsAndSitesMsg(rLeftMolName,
			 rLeftSiteName,
			 rRightMolName,
			 rRightSiteName));
    }

    static
    noKineticConstsXcpt
    fullMessage(const std::string& rLeftMolName,
		const std::string& rLeftSiteName,
		const std::string& rLeftSiteShapeName,
		const std::string& rRightMolName,
		const std::string& rRightSiteName,
		const std::string& rRightSiteShapeName)
    {
      return
	noKineticConstsXcpt
	(mkFullMsg(rLeftMolName,
		   rLeftSiteName,
		   rLeftSiteShapeName,
		   rRightMolName,
		   rRightSiteName,
		   rRightSiteShapeName));
    }
  };
}

#endif // CPX_NOKINETICCONSTSXCPT_H

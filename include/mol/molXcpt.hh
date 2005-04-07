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

#ifndef MOLXCPT_H
#define MOLXCPT_H

#include "domUtils/domUtils.hh"
#include "mzr/mzrXcpt.hh"


namespace bnd
{
  class bindingSite;
  class mol;
  class modification;
  class modMol;
  
  class duplicateSiteNameXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const mol& rMol,
	  const std::string& rDuplicateName);
  public:
    duplicateSiteNameXcpt(const mol& rMol,
			  const std::string& rDuplicateName) :
      mzr::mzrXcpt(mkMsg(rMol,
			 rDuplicateName))
    {}
  };

  class duplicateModNameXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  std::string badModName);
  public:
    duplicateModNameXcpt(xmlpp::Node* pOffendingNode,
			 std::string badModName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
		    badModName))
    {}
  };

  class duplicateModSiteNameXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const std::string& rBadModSiteName);
  public:
    duplicateModSiteNameXcpt(const std::string& rBadModSiteName) :
      mzr::mzrXcpt(mkMsg(rBadModSiteName))
    {}
  };

  class unknownModXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& badModName);

    static std::string
    mkMsgNoElt(const std::string& rBadModName);

  public:
    unknownModXcpt(xmlpp::Node* pOffendingNode,
		   const std::string& rBadModName) throw() :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
		    rBadModName))
    {}

    // This version is used in a utility function of molUnit
    // to find a mod by name, not necessarily in the context
    // of DOM parsing.
    unknownModXcpt(const std::string& rBadModName) throw() :
      mzr::mzrXcpt(mkMsgNoElt(rBadModName))
    {}
  };

  class unknownSiteXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rBadBindingSiteName);
  public:
    unknownSiteXcpt(const xmlpp::Node* pOffendingNode,
		    const std::string& rBadBindingSiteName) throw() :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
		    rBadBindingSiteName))
    {}
  };

  class unknownMolXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rBadMolName);
  public:
    unknownMolXcpt(const xmlpp::Node* pOffendingNode,
		   const std::string& rBadMolName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
		    rBadMolName))
    {}
  };

  class unknownSiteShapeXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const bindingSite& rBindingSite,
	  const mol* pMol,
	  const std::string& rBadSiteShapeName);
  public:
    unknownSiteShapeXcpt(const xmlpp::Node* pOffendingNode,
			 const bindingSite& rBindingSite,
			 const mol* pMol,
			 const std::string& rBadSiteShapeName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
		    rBindingSite,
		    pMol,
		    rBadSiteShapeName))
    {}
  };

  class unknownModSiteXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkDomMsg(xmlpp::Node* pOffendingNode,
	     const std::string& rModSiteName,
	     const modMol* pModMol);

    static std::string
    mkPlainMsg(const std::string& rModSiteName);
  public:
    unknownModSiteXcpt(xmlpp::Node* pOffendingNode,
		       const std::string& rModSiteName,
		       const modMol* pModMol) :
      mzr::mzrXcpt(mkDomMsg(pOffendingNode,
			    rModSiteName,
			    pModMol))
    {}

    unknownModSiteXcpt(const std::string& rModSiteName) :
      mzr::mzrXcpt(mkPlainMsg(rModSiteName))
    {}
  };

  class badModMolCastXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const mol* pMol);
  public:
    badModMolCastXcpt(xmlpp::Node* pOffendingNode,
		      const mol* pMol) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 pMol))
    {}
  };

  class badSmallMolCastXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const mol* pMol);
  public:
    badSmallMolCastXcpt(xmlpp::Node* pOffendingNode,
			const mol* pMol) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 pMol))
    {}
  };

  class badMolParamClassXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const std::string& rParamClassName,
	  const std::string& rMolName);
  public:
    badMolParamClassXcpt(const std::string& rParamClassName,
			 const std::string& rMolName) :
      mzr::mzrXcpt(mkMsg(rParamClassName,
			 rMolName))
    {}
  };

  class modMolQueryTargetXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const modification* pMod,
	  int modNdx);
  public:
    modMolQueryTargetXcpt(const modification* pMod,
			  int modNdx) :
      mzr::mzrXcpt(mkMsg(pMod,
			 modNdx))
    {}
  };
}

#endif // MOLXCPT_H

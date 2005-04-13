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
#include "mol/molXcpt.hh"
#include "mol/mol.hh"
#include "mol/modMol.hh"

namespace bnd
{
  std::string
  duplicateSiteNameXcpt::mkMsg(const mol& rMol,
			       const std::string& rDuplicateName)
  {
    std::ostringstream msgStream;
    msgStream << "Attempted to add site with duplicate name `"
	      << rDuplicateName
	      << "' to mol "
	      << rMol.getName();
    return msgStream.str();
  }

  std::string
  duplicateModNameXcpt::mkMsg(xmlpp::Node* pOffendingNode,
			      std::string badModName)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Name of modification `"
	      << badModName
	      << "' duplicates another modification name.";
    return msgStream.str();
  }

  std::string
  duplicateModSiteNameXcpt::mkMsg(const std::string& rBadModSiteName)
  {
    std::ostringstream msgStream;
    msgStream << "Duplicate modification site name `"
	      << rBadModSiteName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownModXcpt::mkMsg(xmlpp::Node* pOffendingNode,
			const std::string& badModName)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Unknown modification `"
	      << badModName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownModXcpt::mkMsgNoElt(const std::string& rBadModName)
  {
    std::ostringstream msgStream;
    msgStream << "Unknown modification `"
	      << rBadModName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownSiteXcpt::mkMsg(const xmlpp::Node* pOffendingNode,
			 const std::string& rBadBindingSiteName)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Unknown binding site `"
	      << rBadBindingSiteName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownMolXcpt::mkMsg(const xmlpp::Node* pOffendingNode,
			const std::string& rBadMolName)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Unknown mol `"
	      << rBadMolName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownSiteShapeXcpt::mkMsg(const xmlpp::Node* pOffendingNode,
			      const bindingSite& rBindingSite,
			      const mol* pMol,
			      const std::string& rBadSiteShapeName)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Binding site "
	      << rBindingSite.getName()
	      << " on mol "
	      << pMol->getName()
	      << " has no shape named `"
	      << rBadSiteShapeName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownModSiteXcpt::mkDomMsg(xmlpp::Node* pOffendingNode,
			       const std::string& rModSiteName,
			       const modMol* pModMol)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Mod-mol "
	      << pModMol->getName()
	      << " has no modification site named `"
	      << rModSiteName
	      << "'.";
    return msgStream.str();
  }

  std::string
  unknownModSiteXcpt::mkPlainMsg(const std::string& rModSiteName)
  {
    std::ostringstream msgStream;
    msgStream << "Mod-mol (name unknown) "
	      << "has no modification site named `"
	      << rModSiteName
	      << "'.";;
    return msgStream.str();
  }

  std::string
  badModMolCastXcpt::mkMsg(xmlpp::Node* pOffendingNode,
			   const mol* pMol)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Mol "
	      << pMol->getName()
	      << " is not a mod-mol.";
    return msgStream.str();
  }

  std::string
  badSmallMolCastXcpt::mkMsg(xmlpp::Node* pOffendingNode,
			     const mol* pMol)
  {
    std::ostringstream msgStream;
    msgStream << domUtils::domXcptMsg(pOffendingNode)
	      << "Mol "
	      << pMol->getName()
	      << " is not a small-mol.";
    return msgStream.str();
  }

  std::string
  badMolParamClassXcpt::mkMsg(const std::string& rParamClassName,
			      const std::string& rMolName)
  {
    std::ostringstream msgStream;
    msgStream << mzr::internalXcptMsg()
	      << "Parameter was of the wrong class, "
	      << rParamClassName
	      << ", for mol "
	      << rMolName
	      << ".";
    return msgStream.str();
  }

  std::string
  modMolQueryTargetXcpt::mkMsg(const modification* pMod,
			       int modNdx)
  {
    std::ostringstream msgStream;
    msgStream << mzr::internalXcptMsg()
	      << "Could not test for modification "
	      << pMod->name
	      << " at index "
	      << modNdx
	      << " since mol is not a mod-mol.";
    return msgStream.str();
  }
}

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

#ifndef PLEXXCPT_H
#define PLEXXCPT_H

#include "domUtils/domUtils.hh"
#include "mzr/mzrXcpt.hh"

namespace plx
{
  class plexFamily;
  class parserPlex;

  // Need to change all of these to emit standard line number, etc.

  class unknownMolInstanceXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rMolInstName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Unknown mol instance `"
		<< rMolInstName
		<< ".";
      return msgStream.str();
    }
  public:
    unknownMolInstanceXcpt(const xmlpp::Node* pOffendingNode,
			   const std::string& rMolInstName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 rMolInstName))
    {}
  };

  // Maybe it would be better to have a general exception that would
  // be useful for all mol types?
  class instanceNotModMolXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rInstanceName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Mol instance "
		<< rInstanceName
		<< " is not a mod-mol.";
      return msgStream.str();
    }
  public:
    instanceNotModMolXcpt(const xmlpp::Node* pOffendingNode,
			  const std::string& rInstanceName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 rInstanceName))
    {}
  };

  class multiplyBoundSiteXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rBindingSiteName,
	  const std::string& rMolInstanceName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Binding site "
		<< rBindingSiteName
		<< " on mol instance "
		<< rMolInstanceName
		<< " is in more than one binding.";
      return msgStream.str();
    }
  public:
    multiplyBoundSiteXcpt(const xmlpp::Node* pOffendingNode,
			  const std::string& rBindingSiteName,
			  const std::string& rMolInstanceName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 rBindingSiteName,
			 rMolInstanceName))
    {}
  };

  // Thrown in plexAllostery.cc, plexFamily::insertFreeSiteParam::operator().
  class badFreeSiteParamXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "Could not find free site parameter in allostery map.";
      return msgStream.str();
    }
  public:
    badFreeSiteParamXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {}
  };

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
  class noKineticConstsXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const std::string& rLeftSiteShapeName,
	  const std::string& rRightSiteShapeName)
    {
      std::ostringstream msgStream;
      msgStream << "Attempt to form binding of site shapes "
		<< rLeftSiteShapeName
		<< " and "
		<< rRightSiteShapeName
		<< " for which no on/off rates have been given.";
      return msgStream.str();
    }

    static std::string
    mkFeatureMsg(const std::string& rLeftMolName,
		 const std::string& rLeftSiteName,
		 const std::string& rRightMolName,
		 const std::string& rRightSiteName)
    {
      std::ostringstream msgStream;
      msgStream << "Plex contains binding between site "
		<< rLeftSiteName
		<< " on mol "
		<< rLeftMolName
		<< " and site "
		<< rRightSiteName
		<< " on mol "
		<< rRightMolName
		<< " for which no on/off rates have been given.";
      return msgStream.str();
    }

    // Thrown in dimerizeRxnGen.cc, makeBinaryReactions.
    static std::string
    mkBindingMsg(const std::string& rLeftMolName,
		 const std::string& rLeftSiteName,
		 const std::string& rLeftSiteShapeName,
		 const std::string& rRightMolName,
		 const std::string& rRightSiteName,
		 const std::string& rRightSiteShapeName)
    {
      std::ostringstream msgStream;
      msgStream << "Attempt to generate binding reaction between site "
		<< rLeftSiteName
		<< " on mol "
		<< rLeftMolName
		<< " (in shape "
		<< rLeftSiteShapeName
		<< ") with site "
		<< rRightSiteName
		<< " on mol "
		<< rRightMolName
		<< " (in shape "
		<< rRightSiteShapeName
		<< ") for which no on/off rates have been given.";
      return msgStream.str();
    }
      
  public:
    noKineticConstsXcpt(const std::string& rLeftSiteShapeName,
			const std::string& rRightSiteShapeName) :
      mzr::mzrXcpt(mkMsg(rLeftSiteShapeName,
			 rRightSiteShapeName))
      {}

    noKineticConstsXcpt(const std::string& rLeftMolName,
			const std::string& rLeftSiteName,
			const std::string& rRightMolName,
			const std::string& rRightSiteName) :
      mzr::mzrXcpt(mkFeatureMsg(rLeftMolName,
				rLeftSiteName,
				rRightMolName,
				rRightSiteName))
    {}

    noKineticConstsXcpt(const std::string& rLeftMolName,
			const std::string& rLeftSiteName,
			const std::string& rLeftSiteShapeName,
			const std::string& rRightMolName,
			const std::string& rRightSiteName,
			const std::string& rRightSiteShapeName) :
      mzr::mzrXcpt(mkBindingMsg(rLeftMolName,
				rLeftSiteName,
				rLeftSiteShapeName,
				rRightMolName,
				rRightSiteName,
				rRightSiteShapeName))
    {}
  };

  // Thrown in plexAllostery.cc, plexFamily::insertBindingParam::operator().
  class missingLeftParamForBindingXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "left site missing from plex allostery database.";
      return msgStream.str();
    }
  public:
    missingLeftParamForBindingXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {}
  };

  // Thrown in plexAllostery.cc, plexFamily::insertBindingParam::operator().
  class missingRightParamForBindingXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "right site missing from plex allostery database.";
      return msgStream.str();
    }
  public:
    missingRightParamForBindingXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {}
  };

  class incorrectFlipFlagsXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<<"Got incorrect flip flag while constructing plexIsoPair.";
      return msgStream.str();
    }
  public:
    incorrectFlipFlagsXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {}
  };

  class plexNotConnectedXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<<"Left plex not connected in plex isomorphism search.";
      return msgStream.str();
    }
  public:
    plexNotConnectedXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {}
  };

  // Thrown in cxSiteParam.hh, cxSite.getSiteParam().
  class badSiteSpecXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "Invalid site spec for complex.";
      return msgStream.str();
    }
  public:
    badSiteSpecXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {
    }
  };

  // Thrown in siteToShapeMap::getSiteShape.
  class unmappedSiteSpecXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "There is no binding site with the given specification.";
      return msgStream.str();
    }
  public:
    unmappedSiteSpecXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {
    }
  };

  // Thrown in plexFamily::mustFindOmniForNode.
  class noOmniForNodeXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pParentNode)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< domUtils::domXcptMsg(pParentNode)
		<< "No omniPlex was found for this node. Unregistered "
		<< "omniPlex Xpath?";
      return msgStream.str();
    }
  public:
    noOmniForNodeXcpt(xmlpp::Node* pParentNode) :
      mzr::mzrXcpt(mkMsg(pParentNode))
    {}
  };
}

#endif // PLEXXCPT_H

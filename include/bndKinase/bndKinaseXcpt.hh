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

#ifndef BNDKINASEXCPT_H
#define BNDKINASEXCPT_H

#include "domUtils/domUtils.hh"
#include "mzr/mzrXcpt.hh"

namespace bndKinase
{
  class badSmallMolInstanceXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rInstanceName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "expected mol instance "
		<< rInstanceName
		<< " to be a small-mol, but it was not.";
      return msgStream.str();
    }
  public:
    badSmallMolInstanceXcpt(const xmlpp::Node* pOffendingNode,
		    const std::string& rInstanceName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 rInstanceName))
    {}
  };

  class badSmallMolXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rMolName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "expected mol "
		<< rMolName
		<< " to be a small-mol, but it was not.";
      return msgStream.str();
    }
  public:
    badSmallMolXcpt(const xmlpp::Node* pOffendingNode,
		    const std::string& rMolName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 rMolName))
    {}
  };

  class badModMolInstanceXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pOffendingNode,
	  const std::string& rMolInstanceName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "expected mol instance "
		<< rMolInstanceName
		<< " to be a mod-mol, but it was not.";
      return msgStream.str();
    }
  public:
    badModMolInstanceXcpt(const xmlpp::Node* pOffendingNode,
			  const std::string rMolInstanceName) :
      mzr::mzrXcpt(mkMsg(pOffendingNode,
			 rMolInstanceName))
    {}
  };

  class exchangedModXcpt : public mzr::mzrXcpt
  {
    std::string
    mkMsg(bnd::mol* pOffendingMol)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "Expected mod mol, but got mol "
		<< pOffendingMol->getName()
		<< ".";
      return msgStream.str();
    }
  public:
    exchangedModXcpt(bnd::mol* pOffendingMol) :
      mzr::mzrXcpt(mkMsg(pOffendingMol))
    {}
  };
}

#endif // BNDKINASEXCPT_H

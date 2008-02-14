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
#include "utl/badChildCountXcpt.hh"

namespace utl
{
  namespace dom
  {
    std::string
    badChildCountXcpt::
    mkGeneralMsg(const xmlpp::Node* pParentNode,
		 const std::string& rChildName,
		 int requiredCount,
		 int actualCount)
    {
      std::ostringstream msgStream;
      msgStream << xcpt::mkMsg(pParentNode)
		<< "Expected "
		<< requiredCount
		<< " "
		<< rChildName
		<< " elements in "
		<< pParentNode->get_name()
		<< " element; got "
		<< actualCount
		<< ".";
      return msgStream.str();
    }

    std::string
    badChildCountXcpt::
    mkChoiceMsg(const xmlpp::Node* pParentNode,
		int actualCount)
    {
      std::ostringstream msgStream;
      msgStream << xcpt::mkMsg(pParentNode)
		<< "Expected one child of "
		<< pParentNode->get_name()
		<< " node among alternatives "
		<< "offered by choice, got "
		<< actualCount
		<< ".";
      return msgStream.str();
    }

    std::string
    badChildCountXcpt::
    mkOneOrMoreMsg(const xmlpp::Node* pParentNode)
    {
      std::ostringstream msgStream;
      msgStream << xcpt::mkMsg(pParentNode)
		<< "Expected one or more children of "
		<< pParentNode->get_name()
		<< " node among alternatives, got none.";
      return msgStream.str();
    }

    // In support of RNG's "optional" construct.
    std::string
    badChildCountXcpt::
    mkZeroOrOneMsg(const xmlpp::Node* pParentNode,
		   const std::string& rChildName,
		   int actualCount)
    {
      std::ostringstream msgStream;
      msgStream << xcpt::mkMsg(pParentNode)
		<< "Expected zero or one "
		<< rChildName
		<< " elements in "
		<< pParentNode->get_name()
		<< " element; got "
		<< actualCount
		<< ".";
      return msgStream.str();
    }

    badChildCountXcpt::
    badChildCountXcpt(const std::string& rMsg) :
      xcpt(rMsg)
    {}

    badChildCountXcpt
    badChildCountXcpt::
    general(const xmlpp::Node* pParentNode,
	    const std::string& rChildName,
	    int requiredCount,
	    int actualCount)
      throw()
    {
      return badChildCountXcpt(mkGeneralMsg(pParentNode,
					    rChildName,
					    requiredCount,
					    actualCount));
    }

    badChildCountXcpt
    badChildCountXcpt::
    choice(const xmlpp::Node* pParentNode,
	   int actualCount)
      throw()
    {
      return badChildCountXcpt(mkChoiceMsg(pParentNode,
					   actualCount));
    }

    badChildCountXcpt
    badChildCountXcpt::
    oneOrMore(const xmlpp::Node* pParentNode)
      throw()
    {
      return badChildCountXcpt(mkOneOrMoreMsg(pParentNode));
    }

    badChildCountXcpt
    badChildCountXcpt::
    zeroOrOne(const xmlpp::Node* pParentNode,
	      const std::string& rChildName,
	      int actualCount)
      throw()
    {
      return badChildCountXcpt(mkZeroOrOneMsg(pParentNode,
					      rChildName,
					      actualCount));
    }
  }
}

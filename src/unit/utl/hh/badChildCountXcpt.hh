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

#ifndef UTL_BADCHILDCOUNTXCPT_H
#define UTL_BADCHILDCOUNTXCPT_H

#include "utl/dom.hh"

namespace utl
{
  namespace dom
  {
    class badChildCountXcpt : 
      public xcpt
    {
      static std::string
      mkGeneralMsg(const xmlpp::Node* pParentNode,
		   const std::string& rChildName,
		   int requiredCount,
		   int actualCount);

      static std::string
      mkChoiceMsg(const xmlpp::Node* pParentNode,
		  int actualCount);

      static std::string
      mkOneOrMoreMsg(const xmlpp::Node* pParentNode);

      // In support of RNG's "optional" construct.
      static std::string
      mkZeroOrOneMsg(const xmlpp::Node* pParentNode,
		     const std::string& rChildName,
		     int actualCount);

      // This private constructor arranges for creating messages
      // in different ways, returning the same class of exception,
      // under different circumstances.
      badChildCountXcpt(const std::string& rMsg);

    public:

      // For when a definite number of children with a particular name
      // (e.g. 1) is required, but another number appears.
      static
      badChildCountXcpt
      general(const xmlpp::Node* pParentNode,
	      const std::string& rChildName,
	      int requiredCount,
	      int actualCount)
	throw();

      // When there are several possibilities for the child element's name, but
      // only one must appear, as in an RNG schema "choice" construct.
      static
      badChildCountXcpt
      choice(const xmlpp::Node* pParentNode,
	     int actualCount)
	throw();

      // When there are several possibilities for the child element's name, and
      // at least one must appear, as in an RNG schema "oneOrMore" construct.
      static
      badChildCountXcpt
      oneOrMore(const xmlpp::Node* pParentNode)
	throw();

      // Using differently named static functions avoids the problem of
      // constructors all needing different signatures.  When it's time to use
      // it, this looks like throw badChildCountXcpt::optional(...).
      static
      badChildCountXcpt
      zeroOrOne(const xmlpp::Node* pParentNode,
		const std::string& rChildName,
		int actualCount)
	throw();
    };
  }
}

#endif // UTL_BADCHILDCOUNTXCPT_H

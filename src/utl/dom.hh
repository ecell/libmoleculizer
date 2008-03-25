 /////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008  Walter Lawrence (Larry) Lok.
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

#ifndef UTL_DOM_H
#define UTL_DOM_H

#include <string>
#include <set>
#include <libxml++/libxml++.h>
#include "utl/xcpt.hh"

namespace utl
{
  namespace dom
  {
    // Base class for exceptions thrown by utl::dom routines.
    class xcpt :
      public utl::xcpt
    {
    public:
      xcpt(const std::string& rMessage) throw() :
	utl::xcpt(rMessage)
      {}

      xcpt(const char* pMessage) throw() :
	utl::xcpt(pMessage)
      {}

      // Generates "base" diagnostic string for dom error messages, giving the
      // line number and xpath to the the "offending" node.
      //
      // If the default null Node pointer is given, an empty string is
      // returned.
      static std::string
      mkMsg(const xmlpp::Node* pOffendingNode = 0);
    };

    xmlpp::Element*
    mustBeElementPtr(xmlpp::Node* pNode)
      throw(xcpt);

    const xmlpp::Element*
    mustBeElementPtr(const xmlpp::Node* pNode)
      throw(xcpt);

    xmlpp::Element*
    mustGetUniqueChild(const xmlpp::Node* pParentNode,
		       const std::string& rChildName)
      throw(xcpt);

    xmlpp::Element*
    mustGetUniqueChild(const xmlpp::Node* pParentNode,
		       const std::set<std::string>& rChoiceEltNames)
      throw(xcpt);

    // Returns null pointer if no child has the given name.  Throws
    // badChildCountXcpt if more than one child has the given name.
    xmlpp::Element*
    getOptionalChild(const xmlpp::Node* pParentNode,
		     const std::string& rChildName)
      throw(xcpt);

    std::string
    mustGetAttrString(const xmlpp::Element* pElement,
		      const std::string& rAttrName)
      throw(xcpt);

    double
    mustGetAttrDouble(const xmlpp::Element* pElement,
		      const std::string& rAttrName)
      throw(xcpt);

    double
    mustGetAttrPosDouble(const xmlpp::Element* pElement,
			 const std::string& AttrName)
      throw(xcpt);

    double
    mustGetAttrNNDouble(const xmlpp::Element* pElement,
			const std::string& rAttrName)
      throw(xcpt);

    int
    mustGetAttrInt(const xmlpp::Element* pElement,
		   const std::string& rAttrName)
      throw(xcpt);

    int
    mustGetAttrPosInt(const xmlpp::Element* pElement,
		      const std::string& rAttrName)
      throw(xcpt);

    int
    mustGetAttrNNInt(const xmlpp::Element* pElement,
		     const std::string& rAttrName)
      throw(xcpt);

    // Inserts a stereotyped double-valued parameter element, with additional
    // elements giving the parameter value in scientific notation
    //
    // This is basically part of a fix to introduce scientific notation for
    // all parameter values so as to be able, using XSLT, to generate SBML and
    // other formats that have the fraction and exponent in separate XML
    // constructs.
    void
    addDoubleParamChild(xmlpp::Node* pParentNode,
			const std::string& rChildName,
			const std::string& rParameterName,
			double parameterValue);
  }
}

#endif // UTL_DOM_H

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

#ifndef DOMUTILS_H
#define DOMUTILS_H

#include <string>
#include <sstream>
#include <set>
#include <libxml++/libxml++.h>
#include "mzr/mzrXcpt.hh"

namespace domUtils
{
  // Utility function to give line number and xpath to offending
  // element; a sort of "base class mkMsg" for exceptions in dom
  // parsing context (where node is available.)
  std::string
  domXcptMsg(const xmlpp::Node* pOffendingNode);

  // Used in domBatchJob.hh
  class noDocumentParsedXcpt : public mzr::mzrXcpt
  {
  public:
    noDocumentParsedXcpt(void) :
      mzrXcpt("Test of parser shows no document has been parsed.")
    {}
  };

  // Used in domBatchJob.hh
  class insufficientArgsXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(int actualArgCount,
	  int minimumArgCount)
    {
      std::ostringstream msgStream;
      msgStream << "Expected "
		<< minimumArgCount
		<< " command line arguments; got "
		<< actualArgCount
		<< ".";
      return msgStream.str();
    }
  public:
    insufficientArgsXcpt(int actualArgCount,
			 int minimumArgCount) :
      mzr::mzrXcpt(mkMsg(actualArgCount,
			 minimumArgCount))
    {}
  };

  class badElementCastXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pNode)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pNode)
		<< "Could not cast "
		<< pNode->get_name()
		<< " node to xmlpp::Element.";
      return msgStream.str();
    }
  public:
    badElementCastXcpt(const xmlpp::Node* pNode) throw() :
      mzrXcpt(mkMsg(pNode))
    {}
  };

  class missingAttrXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Element* pElt,
	  const std::string& rAttrName)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pElt)
		<< pElt->get_name()
		<< " element has no "
		<< rAttrName
		<< " attribute.";
      return msgStream.str();
    }
  public:
    missingAttrXcpt(const xmlpp::Element* pElt,
		    const std::string& rAttrName) throw() :
      mzr::mzrXcpt(mkMsg(pElt,
		    rAttrName))
    {}
  };

  class badChildCountXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pParentNode,
	  const std::string& rChildName,
	  int requiredCount,
	  int actualCount)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pParentNode)
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

    static std::string
    mkChoiceMsg(const xmlpp::Node* pParentNode,
		int actualCount)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pParentNode)
		<< "Expected one child of "
		<< pParentNode->get_name()
		<< " node among alternatives "
		<< "offered by choice, got "
		<< actualCount
		<< ".";
      return msgStream.str();
    }

    static std::string
    mkOneOrMoreMsg(const xmlpp::Node* pParentNode)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pParentNode)
		<< "Expected one or more children of "
		<< pParentNode->get_name()
		<< " node among alternatives, got none.";
      return msgStream.str();
    }

    // In support of RNG's "optional" construct.
    static std::string
    mkZeroOrOneMsg(const xmlpp::Node* pParentNode,
		  const std::string& rChildName,
		  int actualCount)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pParentNode)
		<< "Expected zero or one "
		<< rChildName
		<< " elements in "
		<< pParentNode->get_name()
		<< " element; got "
		<< actualCount
		<< ".";
      return msgStream.str();
    }

    // This private constructor arranges for creating messages
    // in different ways, returning the same class of exception,
    // under different circumstances.
    badChildCountXcpt(const std::string& rMsg) :
      mzr::mzrXcpt(rMsg)
    {}

  public:

    // Using a different constructor for each different reason for throwing an
    // exception is clearly untenable.  The better way is as for "zeroOrOne"
    // below.

    // When there is only one possibility for the child element's name.
    badChildCountXcpt(const xmlpp::Node* pParentNode,
		      const std::string& rChildName,
		      int requiredCount,
		      int actualCount) throw() :
      mzr::mzrXcpt(mkMsg(pParentNode,
		    rChildName,
		    requiredCount,
		    actualCount))
    {}

    // When there are several possibilities for the child element's name, but
    // only one must appear, as in an RNG schema "choice" construct.
    badChildCountXcpt(const xmlpp::Node* pParentNode,
		      int actualCount) throw() :
      mzr::mzrXcpt(mkChoiceMsg(pParentNode,
			       actualCount))
    {}

    // When there are several possibilities for the child element's name, and
    // at least one must appear, as in an RNG schema "oneOrMore" construct.
    badChildCountXcpt(const xmlpp::Node* pParentNode) throw() :
      mzr::mzrXcpt(mkOneOrMoreMsg(pParentNode))
    {}

    // Using differently named static functions avoids the problem of
    // constructors all needing different signatures.  When it's time to use
    // it, this looks like throw badChildCountXcpt::optional(...).
    static
    badChildCountXcpt
    zeroOrOne(const xmlpp::Node* pParentNode,
	      const std::string& rChildName,
	      int actualCount)
      throw()
    {
      return badChildCountXcpt(mkZeroOrOneMsg(pParentNode,
					      rChildName,
					      actualCount));
    }
  };

  class badDoubleAttrXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Element* pElement,
	  const std::string& rAttrName,
	  const std::string& rBadAttrValue)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pElement)
		<< "Expected double for value of "
		<< rAttrName
		<< " attribute; got `"
		<< rBadAttrValue
		<< "'.";
      return msgStream.str();
    }
  public:
    badDoubleAttrXcpt(const xmlpp::Element* pElement,
		      const std::string& rAttrName,
		      const std::string& rBadAttrValue) throw() :
      mzr::mzrXcpt(mkMsg(pElement,
		    rAttrName,
		    rBadAttrValue))
    {}
  };

  class badPosDoubleAttrXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Element* pOffendingElement,
	  const std::string& rAttrName,
	  double badDoubleValue)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pOffendingElement)
		<< "Expected "
		<< rAttrName
		<< " attribute to be positive; got "
		<< badDoubleValue
		<< ".";
      return msgStream.str();
    }
  public:
    badPosDoubleAttrXcpt(const xmlpp::Element* pOffendingElement,
			 const std::string& rAttrName,
			 double badDoubleValue) throw() :
      mzr::mzrXcpt(mkMsg(pOffendingElement,
		    rAttrName,
		    badDoubleValue))
    {}
  };

  class badIntAttrXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Element* pOffendingElement,
	  const std::string& rAttrName,
	  const std::string& rBadAttrValue)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pOffendingElement)
		<< "Expected integer for value of "
		<< rAttrName
		<< " attribute; got `"
		<< rBadAttrValue
		<< "'.";
      return msgStream.str();
    }
  public:
    badIntAttrXcpt(const xmlpp::Element* pOffendingElement,
		   const std::string& rAttrName,
		   const std::string& rBadAttrValue) :
      mzr::mzrXcpt(mkMsg(pOffendingElement,
		    rAttrName,
		    rBadAttrValue))
    {}
  };

  class badPosIntAttrXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const xmlpp::Element* pOffendingElt,
	  const std::string& rAttrName,
	  int badAttrValue)
    {
      std::ostringstream msgStream;
      msgStream << domXcptMsg(pOffendingElt)
		<< "Expected positive integer for value of "
		<< rAttrName
		<< " attribute; got "
		<< badAttrValue
		<< ".";
      return msgStream.str();
    }
  public:
    badPosIntAttrXcpt(const xmlpp::Element* pOffendingElt,
		      const std::string& rAttrName,
		      int badAttrValue) :
      mzr::mzrXcpt(mkMsg(pOffendingElt,
		    rAttrName,
		    badAttrValue))
    {}
  };

  // There may be a massively more efficient way to do this,
  // but for now....
  template<class writeableType>
  std::string
  stringify(const writeableType& rThingToStringify)
  {
    std::ostringstream oss;
    oss << rThingToStringify;
    return oss.str();
  }

  bool
  stringIsInt(const std::string& rString,
	      int& rInt);

  bool
  stringIsDouble(const std::string& rString,
		 double& rDouble);
  
  xmlpp::Element*
  mustBeElementPtr(xmlpp::Node* pNode)
    throw(badElementCastXcpt);

  const xmlpp::Element*
  mustBeElementPtr(const xmlpp::Node* pNode)
    throw(badElementCastXcpt);

  xmlpp::Element*
  mustGetUniqueChild(const xmlpp::Node* pParentNode,
		     const std::string& rChildName)
    throw(badChildCountXcpt,
	  badElementCastXcpt);

  xmlpp::Element*
  mustGetUniqueChild(const xmlpp::Node* pParentNode,
		     const std::set<std::string>& rChoiceEltNames)
    throw(badChildCountXcpt,
	  badElementCastXcpt);

  // Returns null pointer if no child has the given name.  Throws
  // badChildCountXcpt if more than one child has the given name.
  xmlpp::Element*
  getOptionalChild(const xmlpp::Node* pParentNode,
		   const std::string& rChildName)
    throw(badChildCountXcpt,
	  badElementCastXcpt);

  std::string
  mustGetAttrString(const xmlpp::Element* pElement,
		    const std::string& rAttrName)
    throw(missingAttrXcpt);

  double
  mustGetAttrDouble(const xmlpp::Element* pElement,
		    const std::string& rAttrName)
    throw(missingAttrXcpt,
	  badDoubleAttrXcpt);

  // Tests for strictly positive.
  double
  mustGetAttrPosDouble(const xmlpp::Element* pElement,
		       const std::string& AttrName)
    throw(missingAttrXcpt,
	  badDoubleAttrXcpt,
	  badPosDoubleAttrXcpt);

  int
  mustGetAttrInt(const xmlpp::Element* pElement,
		 const std::string& rAttrName)
    throw(missingAttrXcpt,
	  badIntAttrXcpt);

  // Tests for strictly positive.
  int
  mustGetAttrPosInt(const xmlpp::Element* pElement,
		    const std::string& rAttrName)
    throw(missingAttrXcpt,
	  badIntAttrXcpt,
	  badPosIntAttrXcpt);

  // Utilities connected with generating, rather than parsing, documents.

  // Tests whether there is a unique child with the specified name.  If there
  // is no such child, then one is inserted and a pointer to it is
  // returned. If there is exactly one such child, a pointer to it is
  // returned, unless it is not an element node (a possibility about which the
  // xmlpp documentation isn't exactly forthcoming) in which case a
  // badElementCastXcpt is thrown.  If more than one child with the specified
  // name already exists, a badChildCountXcpt is thrown.
  xmlpp::Element*
  ensureUniqueChild(xmlpp::Element* pParentElt,
		    const std::string& rChildName)
    throw(badElementCastXcpt,
	  badChildCountXcpt);
}

#endif // DOMUTILS_H

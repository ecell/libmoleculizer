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

#include "domUtils/domUtils.hh"

namespace domUtils
{
  std::string
  domXcptMsg(const xmlpp::Node* pOffendingNode)
  {
    std::ostringstream msgStream;
    msgStream << pOffendingNode->get_line()
	      << ":"
	      << pOffendingNode->get_path()
	      << ": ";
    return msgStream.str();
  }

  bool
  stringIsInt(const std::string& rString,
	      int& rInt)
  {
    const char* start = rString.c_str();
    char* pEnd;
    // Setting the base to 0 here means that strings with C-style base
    // indicators (e.g. 0xFFF) can be read successfully.
    rInt = strtol(start, &pEnd, 0);
    return 0 == *pEnd;
  }

  bool
  stringIsDouble(const std::string& rString,
		 double& rDouble)
  {
    const char* start = rString.c_str();
    char* pEnd;
    rDouble = strtod(start, &pEnd);
    return 0 == *pEnd;
  }

  xmlpp::Element*
  mustBeElementPtr(xmlpp::Node* pNode)
    throw(badElementCastXcpt)
  {
    // The dynamic-ness of this cast is intended for initial debugging.
    //
    // I suspect that the xmlpp::Node::NodeList returned by
    // xmlpp::Element:get_children is always a list of element pointers.
    xmlpp::Element* pElt = dynamic_cast<xmlpp::Element*>(pNode);

    if(0 == pElt) throw badElementCastXcpt(pNode);

    return pElt;
  }

  const xmlpp::Element*
  mustBeElementPtr(const xmlpp::Node* pNode)
    throw(badElementCastXcpt)
  {
    // The dynamic-ness of this cast is intended for initial debugging.
    //
    // I suspect that the xmlpp::Node::NodeList returned by
    // xmlpp::Element:get_children is always a list of element pointers.
    const xmlpp::Element* pElt
      = dynamic_cast<const xmlpp::Element*>(pNode);

    if(0 == pElt) throw badElementCastXcpt(pNode);

    return pElt;
  }

  xmlpp::Element*
  mustGetUniqueChild(const xmlpp::Node* pParentNode,
		     const std::string& rChildName)
    throw(badChildCountXcpt,
	  badElementCastXcpt)
  {
    const xmlpp::Node::NodeList children
      = pParentNode->get_children(rChildName);

    if(1 != children.size())
      throw badChildCountXcpt(pParentNode,
			      rChildName,
			      1,
			      children.size());

    xmlpp::Node* pChildNode = children.front();

    return mustBeElementPtr(pChildNode);
  };

  class insertMatchingChildren :
    public std::unary_function<std::string, void>
  {
    const xmlpp::Node* pParent;
    xmlpp::Node::NodeList& rChildren;
  public:
    insertMatchingChildren(const xmlpp::Node* pParentNode,
			   xmlpp::Node::NodeList& rMatchingChildren) :
      pParent(pParentNode),
      rChildren(rMatchingChildren)
    {
    }

    void
    operator()(const std::string& rChildName) const throw(std::exception)
    {
      xmlpp::Node::NodeList childNodes
	= pParent->get_children(rChildName);

      rChildren.insert(rChildren.end(),
		       childNodes.begin(),
		       childNodes.end());
    }
  };

  // This is in support of RNG's "oneOrMore" construct.
  xmlpp::Node::NodeList
  mustGetSomeChildren(const xmlpp::Node* pParentNode,
		      const std::set<std::string>& rChoiceEltNames)
    throw(badChildCountXcpt)
  {
    xmlpp::Node::NodeList result;

    std::for_each(rChoiceEltNames.begin(),
		  rChoiceEltNames.end(),
		  insertMatchingChildren(pParentNode,
					 result));
    if(result.size() < 1)
      throw badChildCountXcpt(pParentNode);

    return result;
  }

  // This is in support of RNG's "choice" construct.
  xmlpp::Element*
  mustGetUniqueChild(const xmlpp::Node* pParentNode,
		     const std::set<std::string>& rChoiceEltNames)
    throw(badChildCountXcpt,
	  badElementCastXcpt)
  {
    xmlpp::Node::NodeList allMatchingChildren;

    std::for_each(rChoiceEltNames.begin(),
		  rChoiceEltNames.end(),
		  insertMatchingChildren(pParentNode,
					 allMatchingChildren));

    if(allMatchingChildren.size() != 1)
      throw badChildCountXcpt(pParentNode,
			      allMatchingChildren.size());

    return mustBeElementPtr(allMatchingChildren.front());
  }

  // This is in support of RNG's "optional" construct.
  xmlpp::Element*
  getOptionalChild(const xmlpp::Node* pParentNode,
		   const std::string& rChildName)
    throw(badChildCountXcpt,
	  badElementCastXcpt)
  {
    const xmlpp::Node::NodeList children
      = pParentNode->get_children(rChildName);

    // Here trying out a better way to arrange different ways
    // of constructing the same exception for use under different
    // circumstances.
    int childCount = children.size();
    switch(childCount)
      {
      case 0 : return 0;
      case 1 : return mustBeElementPtr(children.front());
      default : throw badChildCountXcpt::zeroOrOne(pParentNode,
						   rChildName,
						   childCount);
      }
  }

  std::string
  mustGetAttrString(const xmlpp::Element* pElement,
		    const std::string& rAttrName)
    throw(missingAttrXcpt)
  {
    xmlpp::Attribute* pAttr
      = pElement->get_attribute(rAttrName);

    if(0 == pAttr) throw missingAttrXcpt(pElement,
					 rAttrName);
    return pAttr->get_value();
  }

  double
  mustGetAttrDouble(const xmlpp::Element* pElement,
		    const std::string& rAttrName)
    throw(missingAttrXcpt, badDoubleAttrXcpt)
  {
    std::string attrString
      = mustGetAttrString(pElement,
			  rAttrName);

    double attrDouble = 0.0;
    if(! stringIsDouble(attrString,
				attrDouble))
      throw badDoubleAttrXcpt(pElement,
			      rAttrName,
			      attrString);

    return attrDouble;
  }

  double
  mustGetAttrPosDouble(const xmlpp::Element* pElement,
		       const std::string& rAttrName)
    throw(missingAttrXcpt,
	  badDoubleAttrXcpt,
	  badPosDoubleAttrXcpt)
  {
    double attrDouble
      = mustGetAttrDouble(pElement,
			  rAttrName);

    if(attrDouble <= 0)
      throw badPosDoubleAttrXcpt(pElement,
				 rAttrName,
				 attrDouble);
    return attrDouble;
  }

  int
  mustGetAttrInt(const xmlpp::Element* pElement,
		 const std::string& rAttrName)
    throw(missingAttrXcpt,
	  badIntAttrXcpt)
  {
    std::string attrString
      = mustGetAttrString(pElement,
			  rAttrName);

    int attrInt = -1;
    if(! stringIsInt(attrString,
		     attrInt))
      throw badIntAttrXcpt(pElement,
			   rAttrName,
			   attrString);
    return attrInt;
  }

  int
  mustGetAttrPosInt(const xmlpp::Element* pElement,
		    const std::string& rAttrName)
    throw(missingAttrXcpt,
	  badIntAttrXcpt,
	  badPosIntAttrXcpt)
  {
    int attrInt
      = mustGetAttrInt(pElement,
		       rAttrName);
    if(attrInt <= 0)
      throw badPosIntAttrXcpt(pElement,
			      rAttrName,
			      attrInt);
    return attrInt;
  }

  xmlpp::Element*
  ensureUniqueChild(xmlpp::Element* pParentElt,
		    const std::string& rChildName)
    throw(badElementCastXcpt,
	  badChildCountXcpt)
  {
    xmlpp::Node::NodeList homonymousChildren
      = pParentElt->get_children(rChildName);

    xmlpp::Element* pChild = 0;
    switch(homonymousChildren.size())
      {
      case 0:
	pChild = pParentElt->add_child(rChildName);
	break;
      case 1:
	pChild = mustBeElementPtr(homonymousChildren.front());
	break;
      default:
	// This will try to generate line number, etc. with results
	// that will probably be somewhere between evil and useless.
	//
	// Possibly I can do some general fix to these exceptions that will
	// avoid the line number when it's not there to be found?
	throw badChildCountXcpt(pParentElt,
				rChildName,
				1,
				homonymousChildren.size());
      }
    return pChild;
  }
}



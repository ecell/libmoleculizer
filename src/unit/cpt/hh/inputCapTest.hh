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

#ifndef CPT_INPUTCAPTEST_H
#define CPT_INPUTCAPTEST_H

#include "utl/dom.hh"
#include "cpt/inputCapabilities.hh"

namespace cpt
{
  class modelNodeNotInCap :
    public std::unary_function<xmlpp::Node*, bool>
  {
    inputCapabilities& rInputCap;
  public:
    modelNodeNotInCap(inputCapabilities& rInputCapabilities) :
      rInputCap(rInputCapabilities)
    {}
    
    bool
    operator()(xmlpp::Node* pNode) const 
      throw(utl::xcpt)
    {
      xmlpp::Element* pElt
	= dynamic_cast<xmlpp::Element*>(pNode);

      // We want to return false if the node is not an element, for example,
      // if the node is a comment.
      //
      // Perhaps other tests might be worthwhile to rule out stray text nodes
      // or whatever else might appear as corrupting matter, but this whole
      // checking process was a little daffy.
      return (pElt && (! rInputCap.handlesModelContentElt(pElt)));
    }
  };

  class explicitSpeciesNodeNotInCap :
    public std::unary_function<xmlpp::Node*, bool>
  {
    inputCapabilities& rInputCap;
  public:
    explicitSpeciesNodeNotInCap(inputCapabilities& rInputCapabilities) :
      rInputCap(rInputCapabilities)
    {}
    
    bool
    operator()(xmlpp::Node* pNode) const 
      throw(utl::xcpt)
    {
      xmlpp::Element* pElt 
	= dynamic_cast<xmlpp::Element*>(pNode);

      // We want to return false if the node is not an element, for example,
      // if the node is a comment.
      return (pElt && (! rInputCap.handlesExplictSpeciesContent(pElt)));
    }
  };

  class speciesStreamNodeNotInCap :
    public std::unary_function<xmlpp::Node*, bool>
  {
    inputCapabilities& rInputCap;
  public:
    speciesStreamNodeNotInCap(inputCapabilities& rInputCapabilities) :
      rInputCap(rInputCapabilities)
    {}
    
    bool
    operator()(xmlpp::Node* pNode) const 
      throw(utl::xcpt)
    {
      xmlpp::Element* pElt 
	= dynamic_cast<xmlpp::Element*>(pNode);

      // We want to return false if the node is not an element, for example,
      // if the node is a comment.
      return (pElt && (! rInputCap.handlesSpeciesStreamsContent(pElt)));
    }
  };

  class eventNodeNotInCap :
    public std::unary_function<xmlpp::Node*, bool>
  {
    inputCapabilities& rInputCap;
  public:
    eventNodeNotInCap(inputCapabilities& rInputCapabilities) :
      rInputCap(rInputCapabilities)
    {}
    
    bool
    operator()(xmlpp::Node* pNode) const 
      throw(utl::xcpt)
    {
      xmlpp::Element* pElt 
	= dynamic_cast<xmlpp::Element*>(pNode);

      // We want to return false if the node is not an element, for example,
      // if the node is a comment.
      return (pElt && (! rInputCap.handlesEventsContent(pElt)));
    }
  };

  class reactionGenNotInCap :
    public std::unary_function<xmlpp::Node*, bool>
  {
    inputCapabilities& rInputCap;
  public:
    reactionGenNotInCap(inputCapabilities& rInputCapabilities) :
      rInputCap(rInputCapabilities)
    {}
    
    bool
    operator()(xmlpp::Node* pNode) const 
      throw(utl::xcpt)
    {
      xmlpp::Element* pElt 
	= dynamic_cast<xmlpp::Element*>(pNode);

      // We want to return false if the node is not an element, for example,
      // if the node is a comment.
      return (pElt && (! rInputCap.handlesReactionGensContent(pElt)));
    }
  };
}

#endif // CPT_INPUTCAPTEST_H

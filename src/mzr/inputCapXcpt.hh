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

#ifndef MZR_INPUTCAPXCPT_H
#define MZR_INPUTCAPXCPT_H

#include "utl/dom.hh"

namespace mzr
{
  class unhandledModelContentXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingModelContentNode)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pOffendingModelContentNode)
		<< "No unit claims to handle this "
		<< pOffendingModelContentNode->get_name()
		<< " node in the model section.";
      return msgStream.str();
    }
  public:
    unhandledModelContentXcpt(xmlpp::Node* pOffendingModelContentNode) :
      utl::xcpt(mkMsg(pOffendingModelContentNode))
    {}
  };

  class unhandledExplicitSpeciesContentXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pBadExplicitSpeciesContentNode)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pBadExplicitSpeciesContentNode)
		<< "No unit claims to handle this "
		<< pBadExplicitSpeciesContentNode->get_name()
		<< " node in the explicit species section.";
      return msgStream.str();
    }
  public:
    unhandledExplicitSpeciesContentXcpt(xmlpp::Node* pBadExplicitSpeciesContentNode) :
      utl::xcpt(mkMsg(pBadExplicitSpeciesContentNode))
    {}
  };
    
  class unhandledSpeciesStreamsContentXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pBadSpeciesStreamsContentNode)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pBadSpeciesStreamsContentNode)
		<< "No unit claims to handle this "
		<< pBadSpeciesStreamsContentNode->get_name()
		<< " node in the species streams section.";
      return msgStream.str();
    }
  public:
    unhandledSpeciesStreamsContentXcpt(xmlpp::Node* pBadSpeciesStreamsContentNode) :
      utl::xcpt(mkMsg(pBadSpeciesStreamsContentNode))
    {}
  };

  class unhandledEventsContentXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingEventsContentNode)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pOffendingEventsContentNode)
		<< "No unit claims to handle this "
		<< pOffendingEventsContentNode->get_name()
		<< " node in the events section.";
      return msgStream.str();
    }

  public:
    unhandledEventsContentXcpt(xmlpp::Node* pOffendingEventsContentNode) :
      utl::xcpt(mkMsg(pOffendingEventsContentNode))
    {}
  };

  class unhandledReactionGenXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pBadReactionGenNode)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pBadReactionGenNode)
		<< "No unit claims to handle this "
		<< pBadReactionGenNode->get_name()
		<< " node in the reaction generators section.";
      return msgStream.str();
    }
  public:
    unhandledReactionGenXcpt(xmlpp::Node* pBadReactionGenNode) :
      utl::xcpt(mkMsg(pBadReactionGenNode))
    {
    }
  };
}

#endif // MZR_INPUTCAPXCPT_H

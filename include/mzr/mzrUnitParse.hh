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

#ifndef MZRUNITPARSE_H
#define MZRUNITPARSE_H

#include <sstream>
#include "domUtils/domUtils.hh"

namespace mzr
{
  class unknownSpeciesXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& rSpeciesName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Unknown species `"
		<< rSpeciesName
		<< "'.";
      return msgStream.str();
    }
  public:
    unknownSpeciesXcpt(xmlpp::Node* pOffendingNode,
		       const std::string& rSpeciesName) :
      mzrXcpt(mkMsg(pOffendingNode,
		    rSpeciesName))
    {}
  };

  class unknownDumpableXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& rDumpableName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Unknown dumpable `"
		<< rDumpableName
		<< "'.";
      return msgStream.str();
    }
  public:
    unknownDumpableXcpt(xmlpp::Node* pOffendingNode,
			const std::string& rDumpableName) :
      mzrXcpt(mkMsg(pOffendingNode,
		    rDumpableName))
    {}
  };

  class unknownDumpStreamXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& rDumpStreamName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Unknown dump-stream `"
		<< rDumpStreamName
		<< "'.";
      return msgStream.str();
    }
  public:
    unknownDumpStreamXcpt(xmlpp::Node* pOffendingNode,
			  const std::string& rDumpStreamName) :
      mzrXcpt(mkMsg(pOffendingNode,
		    rDumpStreamName))
    {}
  };

  class noStopEventWarning : public mzrWarning
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "No stop event has been scheduled.";
      return msgStream.str();
    }
  public:
    noStopEventWarning(xmlpp::Node* pEventsElement) :
      mzrWarning(mkMsg(pEventsElement))
    {}
  };

  class unknownStatStreamXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& rBadStreamName)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Unknown stat-stream `"
		<< rBadStreamName
		<< "'.";
      return msgStream.str();
    }
  public:
    unknownStatStreamXcpt(xmlpp::Node* pOffendingNode,
			  const std::string& rBadStreamName) :
      mzrXcpt(mkMsg(pOffendingNode,
		    rBadStreamName))
    {}
  };
}

#endif // MZRUNITPARSE_H

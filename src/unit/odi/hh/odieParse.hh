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

#ifndef ODIEPARSE_H
#define ODIEPARSE_H

#include "odi/dump.hh"

namespace odie
{
  class duplicateSpeciesEntryXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string duplicateSpeciesName)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pOffendingNode)
		<< "Duplicate entry for species `"
		<< duplicateSpeciesName
		<< "'.";
      return msgStream.str();
    }
  public:
    duplicateSpeciesEntryXcpt(xmlpp::Node* pOffendingNode,
			      const std::string duplicateSpeciesName) :
      utl::xcpt(mkMsg(pOffendingNode,
		      duplicateSpeciesName))
      {}
  };

  class unknownSpeciesXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& rUnknownSpeciesName)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pOffendingNode)
		<< "Unknown species: "
		<< rUnknownSpeciesName
		<< ".";
      return msgStream.str();
    }
  public:
    unknownSpeciesXcpt(xmlpp::Node* pOffendingNode,
		       const std::string& rUnknownSpeciesName) :
      utl::xcpt(mkMsg(pOffendingNode,
		      rUnknownSpeciesName))
    {}
  };

  class unknownDumpableXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string& rUnknownDumpableName)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pOffendingNode)
		<< "Unknown dumpable: "
		<< rUnknownDumpableName
		<< ".";
      return msgStream.str();
    }
    
  public:
    unknownDumpableXcpt(xmlpp::Node* pOffendingNode,
			const std::string& rUnknownDumpableName) :
      utl::xcpt(mkMsg(pOffendingNode,
		      rUnknownDumpableName))
    {}
  };

  class odieSpeciesCatalog :
    public std::map<std::string, int>
  {
  public:
    void
    mustAddUnique(xmlpp::Node* pRequestingNode,
		  const std::string& rSpeciesName,
		  int speciesNdx)
      throw(duplicateSpeciesEntryXcpt)
    {
      if(! insert(std::make_pair(rSpeciesName,
				 speciesNdx)).second)
	throw duplicateSpeciesEntryXcpt(pRequestingNode,
					rSpeciesName);
    }

    int
    mustFindSpecies(xmlpp::Node* pRequestingNode,
		    const std::string& rSpeciesName) const
      throw(unknownSpeciesXcpt)
    {
      const_iterator iEntry = this->find(rSpeciesName);
      if(iEntry == end()) throw unknownSpeciesXcpt(pRequestingNode,
						   rSpeciesName);
      return iEntry->second;
    }
  };

  class odieDumpableCatalog : public std::map<std::string,  odieDumpable>
  {
  public:
    const odieDumpable&
    mustFindDumpable(xmlpp::Node* pRequestingNode,
		     const std::string& rDumpableName) const
    {
      const_iterator iEntry = find(rDumpableName);
      // Maybe don't want to involve libmzr?
      if(iEntry == end()) throw unknownDumpableXcpt(pRequestingNode,
						    rDumpableName);
      return iEntry->second;
    }
  };

  // Holds parsed reaction for later use in construction of
  // differential equation.
  class parserReaction
  {
  public:
    std::map<int, int> substrates;
    std::map<int, int> products;
    double rate;
  };
}

#endif // ODIEPARSE_H

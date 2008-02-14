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

#ifndef DUMP_H
#define DUMP_H

#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <string>
#include "utl/dom.hh"

namespace odie
{
  // Thrown by dumpStream::mustSetStream.
  class badFileNameXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pRequestingNode,
	  const std::string& rBadFileName)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pRequestingNode)
		<< "Could not open file `"
		<< rBadFileName
		<< "'.";
      return msgStream.str();
    }
  public:
    badFileNameXcpt(xmlpp::Node* pRequestingNode,
		    const std::string& rBadFileName) :
      utl::xcpt(mkMsg(pRequestingNode,
		      rBadFileName))
    {}
  };

  class odieDumpable :
    public std::vector<int>
  {
    class totVar;
  public:

    std::string name;

    void
    addVar(int varNdx)
    {
      push_back(varNdx);
    }

    // For dumping one column element.  No tab or newline.
    void
    dump(std::ostream& rOs,
	 const double* pConcentrations) const;

    // For dumping a column header. No tab or newline.
    void
    dumpHeader(std::ostream& rOs) const
    {
      rOs << name;
    }
  };

  class dumpStream :
    public std::vector<odieDumpable>
  {
    class dumpDumpables;
    class dumpDumpableHeader;
  
    std::ostream* pOs;
    std::ofstream* pFileStream;
  
  public:
    dumpStream(void) :
      pOs(0),
      pFileStream(0)
    {}
  
    // Using "-" gets cout; "+" gets cerr; others get file.
    void
    mustSetStream(xmlpp::Node* pRequestingNode,
		  const std::string& fileName) throw(badFileNameXcpt);

    // Open file stream must not be closed in destructor, since
    // dumpStream objects are copied.
    //
    // Instead I have to call this on all the dumpStreams when the
    // run is over.
    void
    close(void)
    {
      if(pFileStream) delete pFileStream;
    }

    void
    addDumpable(const odieDumpable& rDumpable)
    {
      push_back(rDumpable);
    }

    void
    dump(const double* pConcentrations,
	 double curTime) const;

    void
    dumpHeaders(void) const;
  };
}

#endif // DUMP_H    
    
  

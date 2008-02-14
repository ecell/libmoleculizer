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

#include <iostream>
#include "fnd/dumpStream.hh"
#include "fnd/badDumpFileXcpt.hh"

namespace fnd
{
  dumpStream::
  dumpStream(const std::string& rFileName)
    throw(utl::xcpt) :
    fileName(rFileName)
  {
    // Establish what the output stream is, opening the output file
    // if called for.  Note that this file doesn't get closed until
    // the destruction of the dumpStream.
    if(fileName == "-") pOs = &std::cout;
    else if(fileName == "+") pOs = &std::cerr;
    else
      {
	// Open the file for appending, so that if a simulation is continued,
	// the output will also be continued instead of starting over at the
	// second start time.
	pFileStream = new std::ofstream(fileName.c_str(),
					std::ios_base::app);
	if(! (*pFileStream))
	  throw badDumpFileXcpt(fileName);

	pOs = pFileStream;
      }
  }

  void
  dumpStream::
  init(void)
    throw(utl::xcpt)
  {
    std::ostream& rOs = getOstream();

    // The header line is a "comment line" for gnuplot.
    rOs << "#";
    
    iterator iEntry = begin();
    if(iEntry != end())
      {
	basicDmpColumn* pColumn = *iEntry;
	
	pColumn->dumpHeader();

	while(++iEntry != end())
	  {
	    rOs << '\t';
	    pColumn = *iEntry;
	    pColumn->dumpHeader();
	  }
      }
    rOs << std::endl;
  }

  void
  dumpStream::
  doDump(void)
  {
    std::ostream& rOs = getOstream();
    
    iterator iEntry = begin();
    if(iEntry != end())
      {
	basicDmpColumn* pColumn = *iEntry;
	pColumn->doDump();
	
	while(++iEntry != end())
	  {
	    rOs << '\t';
	    pColumn = *iEntry;
	    pColumn->doDump();
	  }
      }
    rOs << std::endl;
  }
}

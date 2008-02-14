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

#include <functional>
#include <algorithm>
#include <iostream>
#include "odi/dump.hh"

namespace odie
{
  class prematureDumpXcpt :
    public utl::xcpt
  {
  public:
    prematureDumpXcpt() :
      utl::xcpt("Attempted dump before stream was opened.")
    {}
  };
      
  class odieDumpable::totVar :
    public std::unary_function<int, void>
  {
    double& tot;
    const double* pConc;

  public:
    totVar(double& rTotal,
	   const double* pConcentrations) :
      tot(rTotal),
      pConc(pConcentrations)
    {}

    void
    operator()(int speciesNdx) const
    {
      tot += pConc[speciesNdx];
    }
  };

  void
  odieDumpable::dump(std::ostream& rOs,
		     const double* pConcentrations) const
  {
    double tot = 0.0;
    for_each(begin(),
	     end(),
	     totVar(tot,
		    pConcentrations));
    rOs << tot;
  }

  void
  dumpStream::mustSetStream(xmlpp::Node* pRequestingNode,
			    const std::string& fileName)
    throw(badFileNameXcpt)
  {
    if(fileName == "-") pOs = &std::cout;
    else if(fileName == "+") pOs = &std::cerr;
    else
      {
	pFileStream = new std::ofstream(fileName.c_str());

	if(! (*pFileStream))
	  throw badFileNameXcpt(pRequestingNode,
				fileName);
	pOs = pFileStream;
      }
  }

  // Open file stream must not be closed in destructor, since
  // dumpStream objects are copied.

  class dumpStream::dumpDumpables :
    public std::unary_function<odieDumpable, void>
  {
    std::ostream& rOs;
    const double* pConc;
  
  public:
    dumpDumpables(std::ostream& rOstream,
		  const double *pConcentrations) :
      rOs(rOstream),
      pConc(pConcentrations)
    {}

    void
    operator()(const odieDumpable& rDumpable) const
    {
      rDumpable.dump(rOs,
		     pConc);
      rOs << '\t';
    }
  };

  void
  dumpStream::dump(const double* pConcentrations,
		   double curTime) const
  {
    if(pOs)
      {
	(*pOs) << curTime
	       << '\t';

	for_each(begin(),
		 end(),
		 dumpDumpables(*pOs,
			       pConcentrations));

	(*pOs) << std::endl;
      }
    else
      throw prematureDumpXcpt();
  }

  class dumpStream::dumpDumpableHeader :
    public std::unary_function<odieDumpable, void>
  {
    std::ostream& rOs;
  public:
    dumpDumpableHeader(std::ostream& rOstream) :
      rOs(rOstream)
    {}

    void
    operator()(const odieDumpable& rDumpable) const
    {
      rDumpable.dumpHeader(rOs);
      rOs << '\t';
    }
  };

  void
  dumpStream::dumpHeaders(void) const
  {
    (*pOs) << "#time\t";

    for_each(begin(),
	     end(),
	     dumpDumpableHeader(*pOs));

    (*pOs) << std::endl;
  }
}

  

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

#include "utl/string.hh"

namespace utl
{
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
}

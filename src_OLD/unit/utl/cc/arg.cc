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
#include "utl/insuffArgsXcpt.hh"
#include "utl/badPosIntArgXcpt.hh"
#include "utl/badNNIntArgXcpt.hh"
#include "utl/badNNDoubleArgXcpt.hh"
#include "utl/arg.hh"

namespace utl
{
  std::string
  mustGetArg(int& rArgc,
	     char**& rArgv)
    throw(xcpt)
  {
    if(rArgc <= 0)
      {
	throw insuffArgsXcpt::general();
      }

    std::string argString(*rArgv);
    ++rArgv;
    --rArgc;

    return argString;
  }
  
  int
  argMustBePosInt(const std::string& rArgString)
    throw(xcpt)
  {
    int argValue = -1;

    if(! (stringIsInt(rArgString,
		      argValue)
	  && (0 < argValue)))
      {
	throw badPosIntArgXcpt(rArgString);
      }

    return argValue;
  }

  int
  argMustBeNNInt(const std::string& rArgString)
    throw(xcpt)
  {
    int argValue = -1;

    if(! (stringIsInt(rArgString,
		      argValue)
	  && (0 <= argValue)))
      {
	throw badNNIntArgXcpt(rArgString);
      }

    return argValue;
  }

  double
  argMustBeNNDouble(const std::string& rArgString)
    throw(xcpt)
  {
    double argValue = -1.0;

    if(! (stringIsDouble(rArgString,
			 argValue)
	  && (0.0 <= argValue)))
      {
	throw badNNDoubleArgXcpt(rArgString);
      }

    return argValue;
  }
}

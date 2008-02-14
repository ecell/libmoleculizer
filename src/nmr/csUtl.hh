/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007  Nathan Addy
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
//   Nathan Addy, Research Associate     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////


#ifndef CSUTL_HH
#define CSUTL_HH

#include <string>
#include <set>
#include <sstream>
#include <iostream>
#include <vector>
#include <utility>

// A couple global functions I have chance to use from time to time.

namespace complexspecies
{

  namespace detail
  {
    
  
    template <class ForwardIter,
	      class OutputIter,
	      class UnaryPred>
    OutputIter copy_if(ForwardIter begin,
		       ForwardIter end,
		       OutputIter dest,
		       UnaryPred f)
    {
      while(begin!=end)
	{
	  if( f(*begin))
	    {
	      *dest=*begin;
	      ++dest;
	    }
	  ++begin;
	}
      return dest;
    }
	  
	    
    template <typename WriteableType>
    std::string
    stringify(const WriteableType& rThingToStringify)
    {
      typename std::ostringstream oss;
      oss << rThingToStringify;
      return oss.str();
    }

  }

}
 


#endif



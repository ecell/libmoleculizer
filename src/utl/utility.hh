/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 Nathan Addy
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
/////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <sstream>

#ifndef UTILITY_HH
#define UTILITY_HH

namespace utl
{
    std::string
    getFileName(int argc,
                char* argv[]);
    
    void tokenize(const std::string& str,
                  std::vector<std::string>& tokens,
                  const std::string& deliminator = " ");
    
    template <class T>
    bool from_string(T& t, 
                     const std::string& s)
    {
        std::istringstream iss(s);
        return !(iss >> t).fail();
    }
    
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

  bool
  stringIsInt(const std::string& rString,
	      int& rInt);

  bool
  stringIsDouble(const std::string& rString,
		 double& rDouble);
  

}

#endif

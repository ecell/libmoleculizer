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


#include "utl/utility.hh"
#include "utl/badFileNameXcpt.hh"

namespace utl
{
  std::string
  getFileName(int argc,
              char* argv[])
  {
    std::string filename;

    for(int i = 0; i != argc; ++i)
      {
        if ( std::string(argv[i]) == std::string("-f") )
          {
            if(i + 1 == argc)
              {
                throw badFileNameXcpt();
              }
            else
              {
                filename = std::string(argv[i + 1]);
                return filename;
              }
          }
      }
    // Nothing was found, so throw an exception.
    throw badFileNameXcpt();
    return 0;
  }


  void tokenize(const std::string& str,
                std::vector<std::string>& tokens,
                const std::string& deliminator)
  {
    const std::string::size_type DELIM_SIZE( deliminator.size() );

    std::string::size_type index(0);

    // If the string starts with a deliminator, ignore it. 
    if ( str.find(deliminator) == 0)
      {
        index += DELIM_SIZE;
      }

    while (true)
      {
        std::string::size_type nextDelimStart = str.find(deliminator, index);

        // Two cases:
        // A deliminator is found
        // A deliminator is not found

        if (nextDelimStart != std::string::npos)
          {
            // everything between index and nextDelimStart should be pushed back
            
            if (index == nextDelimStart)
              {
                tokens.push_back( std::string("") );
              }
            else
              {
                tokens.push_back( std::string(str, index, nextDelimStart - index) );
              }
            index = nextDelimStart + DELIM_SIZE;
          }
        else
          {
            tokens.push_back( std::string(str, index, str.size() - index));
            return;
          }
      }
  }



}

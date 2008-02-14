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

#ifndef UTL_XCPT_H
#define UTL_XCPT_H

#include <iostream>
#include <string>
#include <exception>

namespace utl
{
  class xcpt :
    public std::exception
  {
    std::string message;
  public:
    xcpt(const std::string& rMessage) throw() :
      message(rMessage)
    {}

    xcpt(const char* pMessage) throw() :
      message(pMessage)
    {}

    ~xcpt(void)
      throw()
    {}

    const std::string&
    getMessage(void) const throw()
    {
      return message;
    }

    // Implementation of std::exception::what.
    const char*
    what(void) const throw()
    {
      return message.c_str();
    }
    
    void
    warn(void)
    {
      std::cerr << mkWarnMsg()
		<< getMessage()
		<< std::endl;
    }

    void
    wailAndBail(void)
    {
      std::cerr << mkMsg()
		<< getMessage()
		<< std::endl;
      exit(1);
    }

    // But this doesn't include the application name, and it doesn't
    // look good for getting it in where this message is called for.
    //
    // Most of those cases should be handled with wailAndBail, rather than
    // throw, since they are all essentially debugging loose-ends.
    static std::string
    mkMsg(void)
    {
      return std::string("Internal exception: ");
    }

    static std::string
    mkWarnMsg(void)
    {
      return std::string("Warning: ");
    }
  };
}

#endif // UTL_XCPT_H

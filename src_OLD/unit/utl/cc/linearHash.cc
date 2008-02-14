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

#include "utl/linearHash.hh"

namespace utl
{
  size_t
  linearHash::operator()(const size_t& rData) const
  {
    return (rData * multiplier) + summand;
  }

  // These will need to be adjusted, I expect.  Or maybe not.
  const size_t linearHash::multiplier = 2897564231ul;
  const size_t linearHash::summand = 3248630751ul;

  // This could be somewhat templatized....
  // It also appears to be defined in tauApp.cc.
  class charHashAccum : public std::unary_function<char, void>
  {
    size_t& rValue;
    linearHash lh;
  public:
    charHashAccum(size_t& rHashValue) :
      rValue(rHashValue)
    {
    }

    void
    operator()(char c) const
    {
      rValue = lh(rValue + lh((size_t) c));
    }
  };

  size_t
  linearHash::operator()(const std::string& rString) const
  {
    size_t hashValue;
    std::for_each(rString.begin(),
		  rString.end(),
		  charHashAccum(hashValue));
    return hashValue;
  }
}

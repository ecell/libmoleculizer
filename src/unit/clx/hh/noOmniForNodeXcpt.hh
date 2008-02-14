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

#ifndef CLX_NOOMNIFORNODEXCPT_H
#define CLX_NOOMNIFORNODEXCPT_H

#include "utl/dom.hh"

namespace clx
{
  // Thrown when a client of the clxUnit wants to know what omniplex the
  // plexUnit parsed for one of the client's nodes.  If the clxUnit doesn't
  // know, it probably means that the client unit didn't register an Xpath for
  // the node with the clxUnit, so the clxUnit didn't really parse the
  // omniplex.  (This might be a common error in new unit implementation.)
  class noOmniForNodeXcpt : 
    public utl::xcpt
  {
    static std::string
    mkMsg(const xmlpp::Node* pParentNode = 0);

  public:
    noOmniForNodeXcpt(const xmlpp::Node* pParentNode = 0) :
      utl::xcpt(mkMsg(pParentNode))
    {}
  };
}

#endif // CLX_NOOMNIFORNODEXCPT_H

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

#ifndef UTL_FORBOTH_H
#define UTL_FORBOTH_H

namespace utl
{
  // A binary version of std::for_each.
  template<class inputIter1, class inputIter2, class binaryOp>
  binaryOp
  for_both(inputIter1 first1,
	   inputIter1 last1,
	   inputIter2 first2,
	   binaryOp op)
  {
    while(first1 != last1)
      {
	op(*first1, *first2);
	// Postincrement is costly for some iterators.
	++first1;
	++first2;
      }
    return op;
  }
}

#endif // UTL_FORBOTH_H

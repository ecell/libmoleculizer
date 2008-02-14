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

#ifndef MOL_BADMOLPARAMXCPT_H
#define MOL_BADMOLPARAMXCPT_H

#include "utl/xcpt.hh"

namespace bnd
{
  // Internal exception thrown by alloMols when they are asked to intern a
  // state that is not of the right type for the alloMol.
  //
  // In some sense, this exception can't occur yet, since there is only
  // one kind of alloMol, mzrModMol, in Moleculizer at this time.
  class badMolParamXcpt : 
    public utl::xcpt
  {
    static std::string
    mkMsg(const std::string& rParamClassName,
	  const std::string& rMolName);
  public:
    badMolParamXcpt(const std::string& rParamClassName,
			 const std::string& rMolName) :
      utl::xcpt(mkMsg(rParamClassName,
			 rMolName))
    {}
  };

}

#endif // MOL_BADMOLPARAMXCPT_H

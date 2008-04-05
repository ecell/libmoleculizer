/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008  Nathan Addy
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

#ifndef NMRMANGLERFACTORY_HH
#define NMRMANGLERFACTORY_HH

#include "nameAssembler.hh"
#include "complexSpeciesManglerNames.hh"
#include "nameAssemblers.hh"
#include "noSuchNameManglerXcpt.hh"

namespace nmr
{

  template <typename molT>
  class NameManglerFactory
  {
  public:
    typedef typename nmr::NameAssembler<molT> NameAssemblerType;

    NameAssemblerType*
    create(const std::string& manglerName) throw( NoSuchNameManglerXcpt )
    {
//       if( manglerName == manglernames::basicManglerName )
//         {
//           return new basicNameAssembler<molT>;
//         }
//       else if ( manglerName == manglernames::detailedManglerName )
//         {
//           return new readableNameAssembler<molT>;
//         }
//       else if ( manglerName == manglernames::compactManglerName )
//         {
//           return new MangledNameAssembler<molT>;
//         }
//       else
//         {
//           throw NoSuchNameManglerXcpt( manglerName );
//         }

      if( manglerName == manglernames::compactManglerName )
        {
          return new MangledNameAssembler<molT>;
        }
      else
        {
          throw NoSuchNameManglerXcpt( manglerName );
        }

    }

  };


}


#endif

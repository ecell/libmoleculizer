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

namespace nmr
{

  template <typename molT>
  const NameAssembler<molT>* 
  nmrUnit<molT>::getNameAssembler() const throw( utl::xcpt )
  {
    if (!ptrNameAssembler)
      {
        throw utl::xcpt( "No NameMangler has been set.  Therefore we certainly cannot \"get\" it.");
      }
    else
      {
        return ptrNameAssembler;
      }
  }
    
  template <typename molT>
  void nmrUnit<molT>::setDefaultNameMangler( const std::string& nameManglerName ) throw (NoSuchNameManglerXcpt)
  {
    try
      {
        NameAssembler<MolType>* ptrNA = ptrNameManglerFactory->create(nameManglerName);
        
        if (!ptrNameAssembler)
          {
            throw utl::xcpt( "Unknown error." );
          }
        else
          {
            delete ptrNameAssembler;
            ptrNameAssembler = ptrNA;
          }
      }
    catch( NoSuchNameManglerXcpt xcpt)
      {
        // Catch and release...
        throw xcpt;
      }
  }

}

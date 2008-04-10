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
  nmrUnit<molT>::getNameEncoder() const throw( utl::xcpt )
  {
    if (!ptrNameAssembler)
      {
        throw utl::xcpt( "No NameEncoder has been set.  Therefore we certainly cannot \"get\" it.");
      }
    else
      {
        return ptrNameAssembler;
      }
  }
    
  template <typename molT>
  void nmrUnit<molT>::setDefaultNameEncoder( const std::string& nameEncoderName ) throw (NoSuchNameEncoderXcpt)
  {
    try
      {
        NameAssembler<MolType>* ptrNA = ptrNameEncoderFactory->create(nameEncoderName);
        
        delete ptrNameAssembler;
        ptrNameAssembler = ptrNA;
      }
    catch( NoSuchNameEncoderXcpt xcpt)
      {
        // Catch and release...
        std::cout << "HELP ME HELP ME" << endl;
        throw xcpt;
      }
  }

  template <typename molT>
  void nmrUnit<molT>::parseDomInput(xmlpp::Element* pRootElt,
                                    xmlpp::Element* pModelElt,
                                    xmlpp::Element* pStreamsElt) throw(std::exception)
    {
      try
        {

          // Redo this, as it is much (maybe) worse than it should be.
          std::string namingConventionStyle = utl::dom::mustGetAttrString( pModelElt,
                                                                           eltName::namingConvention);
          setDefaultNameEncoder( namingConventionStyle );

        }
      catch(...)
        {
          return;
        }
    }


  template <typename molT>
  void nmrUnit<molT>::insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
  {
    // For now, do nothing.  It may be best that when the parseDomInput function is 
    // implemented, we will want to insert the generator we are using too.
  }

}

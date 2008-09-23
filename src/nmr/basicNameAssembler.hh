/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008 The Molecular Sciences Institute
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Original Author:
//   Nathan Addy, Research Associate     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//                     
//   
/////////////////////////////////////////////////////////////////////////////

#ifndef __BASICNAMEASSEMBLER_HH
#define __BASICNAMEASSEMBLER_HH

#include "nmrExceptions.hh"
#include "nameAssembler.hh"

#include "utl/macros.hh"

namespace nmr
{

    DECLARE_CLASS( basicNameAssembler);


    // TODO/2 Write a description of the mangling algorithm used by basicNameAssembler::createNameFromOutputState.
    class basicNameAssembler : public NameAssembler
    {
    public:
      std::string createNameFromOutputState( ComplexOutputStateCref aCOS) const;
      ComplexOutputState createOutputStateFromName( const std::string& aMangledName) const throw(utl::NotImplementedXcpt);
    };

}
#endif

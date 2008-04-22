/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008  Nathan Addy
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
//   Nathan Addy, Research Associate     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////


#ifndef __MANGLEDNAMEASSEMBLER_HH
#define __MANGLEDNAMEASSEMBLER_HH

#include "nmrExceptions.hh"
#include "nameAssembler.hh"
#include "complexSpecies.hh"
#include "complexOutputState.hh"
#include "utl/utility.hh"

#include <iostream>
#include <string>
#include <vector>

namespace nmr
{

    class MangledNameAssembler : public NameAssembler
    {
    public:

        typedef std::vector<std::string> strVect;
        typedef std::vector<std::string>::iterator strVectIter;
        typedef std::vector<std::string>::const_iterator cstrVectIter;

        MangledNameAssembler()
            :
            NameAssembler("MangledNameAssembler")
        {
            try
            {
                assertEncodeDecodeAccuracy(this);
            }
            catch(encodeDecodeInconsistencyXcpt xcpt)
            {
                xcpt.warn();
                std::cerr << "Continuing..." << std::endl;
            }
        }

        std::string 
        createNameFromOutputState( const ComplexOutputState& aCOS) const;
        
        ComplexOutputState createOutputStateFromName(const std::string& name) const;

    protected:

        std::string constructMangledMolList(const ComplexOutputState& aComplexOutputState) const;
        std::string constructMangledBindingList(const ComplexOutputState& aComplexOutputState) const;
        std::string constructMangledModificationList(const ComplexOutputState& aComplexOutputState) const;
        std::string getEncodedLength(const std::string& stringInQuestion) const;

        void parseMolString(const std::string& molString, std::vector<ComplexOutputState::MolTokenStr>& molTokenVector) const;
        void parseBindingString(const std::string& bindingString, std::vector<ComplexOutputState::BindingTokenStr>& bindingTokenVector ) const;
        void parseModificationString(const std::string& modString, std::vector<ComplexOutputState::ModificationTokenStr>& modificationTokenVector) const;
   
        std::string processBindingString(const std::string& aString) const
        {
            if (aString.empty())
            {
                throw GeneralNmrXcpt("Error.  MangledNameAssembler::processBindingString failed.\nThis function cannot be called on a null string.");
            }
      

            std::string returnVal;

            if (aString.size()==1)
            {
                returnVal = aString;
            }
            else if (aString.size()<10)
            {
                returnVal = "_" + utl::stringify(aString.size())+ aString;
            }
            else 
            {
                returnVal = "__"+utl::stringify(aString.size()) + "_" + aString;
            }
      
            return returnVal;
        }

        std::string processModificationToken(const std::string& aModString) const;
        

    };
    

}

#endif

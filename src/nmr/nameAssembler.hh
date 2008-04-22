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

#ifndef __NAMEASSEMBLER_HH
#define __NAMEASSEMBLER_HH

#include "complexSpeciesOutputMinimizer.hh"
#include "nmrExceptions.hh"
#include "complexSpecies.hh"
#include "complexOutputState.hh"


#include <iterator>
#include <vector>
#include <string>
#include <utility>

namespace nmr
{

    class NameAssembler
    {
    public:
        NameAssembler(const std::string& name)
            :
            assemblerName(name)
        {}

        virtual ~NameAssembler()
        {}

        const std::string& 
        getName() const 
        {
            return assemblerName;
        }

        std::string 
        createName(ComplexSpeciesCref aComplexSpecies) const
        {
            ComplexOutputState anOutputState;
            aComplexSpecies.constructOutputState( anOutputState );

            return createNameFromOutputState( anOutputState );
        }

        std::string 
        createCanonicalName( ComplexSpeciesCref aComplexSpecies) const
        {
            ComplexSpeciesOutputMinimizer canonicalNameGenerator;
            ComplexOutputState minimalOutputState = canonicalNameGenerator.getMinimalOutputState( aComplexSpecies );
            return createNameFromOutputState( minimalOutputState );
        }

        virtual std::string createNameFromOutputState( ComplexOutputStateCref aComplexOutputState) const = 0;
        virtual ComplexOutputState createOutputStateFromName( const std::string& aComplexOutputState) const = 0;
    
    protected:
        const std::string assemblerName;
        
    protected:
        // This function is intended as a debugish sort of a function that can be called during constructors of 
        // different NameAssemblers.  It will throw an exception if startingComplexOutputState != endingComplexOutputState
        static void assertEncodeDecodeAccuracy(NameAssembler* ptrNameAssembler) throw(encodeDecodeInconsistencyXcpt) 
        {

            // TODO: Create a much better and 2more complicated ComplexOutputState...
            ComplexOutputState aComplexOutputState;

            ComplexOutputState decodedEncodingOS = ptrNameAssembler->createOutputStateFromName( ptrNameAssembler->createNameFromOutputState( aComplexOutputState ) );
            if (! (decodedEncodingOS == aComplexOutputState)  )
            {
                throw encodeDecodeInconsistencyXcpt(ptrNameAssembler->getName());
            }
        }

    };

}

#endif

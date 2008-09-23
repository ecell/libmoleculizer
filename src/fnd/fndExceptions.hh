/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 The Molecular Sciences Institute
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

#ifndef FNDEXCEPTIONS_HH
#define FNDEXCEPTIONS_HH

#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "utl/xcpt.hh"

namespace fnd
{

    class DuplicatedCatalogEntryXcpt : public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& ElementName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: "
                      << "Object " 
                      << "'"
                      << ElementName 
                      << "' was already present in the ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
        DuplicatedCatalogEntryXcpt( const std::string& refObjectName) 
            :
            utl::xcpt( mkMsg( refObjectName ) )
        {}
    };

    class NoSuchSpeciesXcpt : public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& speciesName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: Species '"
                      << speciesName 
                      << "' was not found in ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
        NoSuchSpeciesXcpt( const std::string& speciesName )
            :
            xcpt( mkMsg( speciesName ) )
        {}
        
    };

    class NoSuchReactionXcpt : public utl::xcpt
    {
        std::string
        makeMessage( const std::string& reactionName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: Reaction '"
                      << reactionName 
                      << "' was not found in ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
      ~NoSuchReactionXcpt() throw(){}


        NoSuchReactionXcpt( const std::string& reactionName )
            :
	  xcpt(""),
	  theReactionName( reactionName)
      {
	message = this->makeMessage( reactionName );
      }

    private:
      std::string theReactionName;
    };


}



#endif

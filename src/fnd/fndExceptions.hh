#ifndef FNDEXCEPTIONS_HH
#define FNDEXCEPTIONS_HH

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
        static std::string
        mkMsg( const std::string& reactionName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: Reaction '"
                      << reactionName 
                      << "' was not found in ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
        NoSuchReactionXcpt( const std::string& reactionName )
            :
            xcpt( mkMsg( reactionName ) )
        {}
    };


}



#endif

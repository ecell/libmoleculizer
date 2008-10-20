#include "mzrException.hh"

namespace mzr
{
std::string
dumpableNotSpeciesStreamXcpt::
mkMsg( const std::string& rDumpableName )
{
    std::ostringstream msgStream;
    msgStream << utl::xcpt::mkMsg()
    << "Dumpable '"
    << rDumpableName
    << "' is not a species stream, but "
    << "was registered as such.";
    return msgStream.str();
}

std::string
illegalSpeciesNameXcpt::
mkMsg( const std::string& rSpeciesName,
       const std::string& rMessage )
{
    return ( rSpeciesName + " is an illegal name" );
}

std::string
unkDumpableXcpt::
mkMsg( const std::string& rDumpableName,
       xmlpp::Node* pOffendingNode )
{
    std::ostringstream msgStream;
    msgStream << utl::dom::xcpt::mkMsg( pOffendingNode )
    << "Unknown dumpable `"
    << rDumpableName
    << "'.";
    return msgStream.str();
}

std::string
dupDumpableNameXcpt::
mkMsg( const std::string& rDumpableName,
       xmlpp::Node* pOffendingNode )
{
    std::ostringstream msgStream;
    msgStream << utl::dom::xcpt::mkMsg( pOffendingNode )
    << "There is already a dumpable named "
    << rDumpableName
    << ".";
    return msgStream.str();
}

std::string
unkSpeciesXcpt::
mkMsg( const std::string& rSpeciesName,
       xmlpp::Node* pOffendingNode )
{
    std::ostringstream msgStream;
    msgStream << utl::dom::xcpt::mkMsg( pOffendingNode )
    << "Unknown species `"
    << rSpeciesName
    << "'.";
    return msgStream.str();
}

std::string
dupSpeciesNameXcpt::
mkMsg( const std::string& rSpeciesName,
       xmlpp::Node* pRequestingNode )
{
    std::ostringstream msgStream;
    msgStream << utl::dom::xcpt::mkMsg( pRequestingNode )
    << "There is already a species named "
    << rSpeciesName
    << ".";
    return msgStream.str();
}

std::string
stopEventInPastXcpt::
mkMsg( double now,
       double badEventTime )
{
    std::ostringstream msgStream;
    msgStream << "Stop scheduled at time "
    << badEventTime
    << ", which is in the past at simulation time "
    << now
    << ".";
    return msgStream.str();
}


std::string
unkStatStreamXcpt::
mkMsg( const std::string& rBadStreamName,
       xmlpp::Node* pOffendingNode )
{
    std::ostringstream msgStream;
    msgStream << utl::dom::xcpt::mkMsg( pOffendingNode )
    << "Unknown stat-stream `"
    << rBadStreamName
    << "'.";
    return msgStream.str();
}

std::string
missingExtrapolationParameter::
mkMsg( const std::string& rSpeciesName,
       const std::string& missingParamName )
{
    std::ostringstream msgStream;
    msgStream << "Error: Species '"
    << rSpeciesName
    << "' expected having parameter "
    << missingParamName
    << ", which was not found.";
    return msgStream.str();
}
}

#include "propXcpt.hpp"



std::string
PropertyDoesNotExistXcpt::mkMsg( const std::string& propName)
{
    std::ostringstream oss;
    oss << "Property '" 
        << propName 
        << "' does not exist.";
    return oss.str();
}

std::string
PropertyDoesNotExistXcpt::mkMsg( const std::string& propName,
                                 const PropertiedClass* ptrPropClass)
{
    std::ostringstream oss;
    oss << "Property '" 
        << propName 
        << "' does not exist.";
    return oss.str();
}

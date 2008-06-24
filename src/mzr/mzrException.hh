#ifndef __MZREXCEPTION_HH
#define __MZREXCEPTION_HH

#include "utl/xcpt.hh"

namespace mzr
{

    class BadRulesDefinitionXcpt : public utl::xcpt
    {
    public:
        static std::string 
        mkMsg()
        {
            std::ostringstream oss;
            oss << "Can't parse a bad ruleset.";
            return oss.str();
        }

        BadRulesDefinitionXcpt()
            :
            xcpt(mkMsg())
        {
        }
    };

    
    class IllegalNameXcpt : public utl::xcpt
    {
    public:
        static std::string
        mkMsg( const std::string& illegalName)
        {
            std::ostringstream oss;
            oss << "Illegal Name: '" << illegalName << "'";
            return oss.str();
        }
        
        IllegalNameXcpt( const std::string& illegalName)
            :
            xcpt( mkMsg(illegalName))
        {}
    };

}

#endif

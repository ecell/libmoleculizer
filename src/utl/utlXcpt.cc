// This file contains the implementation mkMsgs for a variety of Larry written 
// exceptions:

// badNNIntArgXcpt
// badPosIntArgXcpt

namespace utl
{
    std::string
    badNNIntArgXcpt::
    mkMsg(const std::string& rTheBadArgument)
    {
        std::ostringstream msgStream;
        msgStream << "Expected non-negative integer for command-line argument; "
        << "got `"
        << rTheBadArgument
        << "'.";
        return msgStream.str();
    }


    std::string
    badPosIntArgXcpt::
    mkMsg(const std::string& rTheBadArgument)
    {
        std::ostringstream msgStream;
        msgStream << "Expected positive integer for command-line argument; got `"
        << rTheBadArgument
        << "'.";
        return msgStream.str();
    }

        std::string
        badPosIntAttrXcpt::
        mkMsg (const xmlpp::Element* pOffendingElt,
               const std::string& rAttrName,
               int badAttrValue)
        {
            std::ostringstream msgStream;
            msgStream << xcpt::mkMsg (pOffendingElt)
            << "Expected positive integer for value of "
            << rAttrName
            << " attribute; got "
            << badAttrValue
            << ".";
            return msgStream.str();
        }

    std::string
    insuffArgsXcpt::
    mkCountsMsg (int actualArgCount,
                 int minimumArgCount)
    {
        std::ostringstream msgStream;
        msgStream << "Expected "
        << minimumArgCount
        << " command line arguments; got "
        << actualArgCount
        << ".";
        return msgStream.str();
    }

    std::string
    insuffArgsXcpt::
    mkGeneralMsg (void)
    {
        return std::string ("Insufficient command-line arguments or "
                            "missing command-line argument.");
    }

    insuffArgsXcpt::
    insuffArgsXcpt (const std::string& rMsg) :
            utl::xcpt (rMsg)
    {}

    insuffArgsXcpt
    insuffArgsXcpt::
    general (void)
    {
        return insuffArgsXcpt (mkGeneralMsg() );
    }

    insuffArgsXcpt
    insuffArgsXcpt::
    counts (int actualArgCount,
            int minimumArgCount)
    {
        return insuffArgsXcpt (mkCountsMsg (actualArgCount,
                                            minimumArgCount) );
    }

}

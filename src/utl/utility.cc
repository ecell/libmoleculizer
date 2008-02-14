#include "utl/utility.hh"
#include "utl/badFileNameXcpt.hh"

namespace utl
{
  std::string
  getFileName(int argc,
              char* argv[])
  {
    std::string filename;

    for(int i = 0; i != argc; ++i)
      {
        if ( std::string(argv[i]) == std::string("-f") )
          {
            if(i + 1 == argc)
              {
                throw badFileNameXcpt();
              }
            else
              {
                filename = std::string(argv[i + 1]);
                return filename;
              }
          }
      }
    // Nothing was found, so throw an exception.
    throw badFileNameXcpt();
    return 0;
  }

}

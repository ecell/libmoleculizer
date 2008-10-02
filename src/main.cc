//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#include "utl/defs.hh"
#include "utl/utlXcpt.hh"
#include "utl/arg.hh"
#include "utl/utility.hh"
#include "utl/linearHash.hh"
#include "mzr/moleculizer.hh"

void
processCommandLineArgs (int argc,
                        char* argv[],
                        mzr::moleculizer& theMzrObject,
                        std::string* fileName);

void displayHelpMessage();


int
main (int argc, char** argv)
{

    try
    {
        std::string filename;
        mzr::moleculizer theApp;

        processCommandLineArgs ( argc, argv, theApp, &filename);
        theApp.attachFileName ( filename );
        
        // Do nothing.

        return 0;
    }
    catch (const utl::xcpt& rException)
    {
        // If a moleculizer exception is thrown report it and crash.
        std::cerr << rException.what() << std::endl;
        return 1;
    }
    catch (const std::exception& rExcept)
    {
        // If an exception is thrown, report it and crash.
        std::cerr << rExcept.what()
        << std::endl;
        return 1;
    }
}

void
processCommandLineArgs (int argc,
                        char* argv[],
                        mzr::moleculizer& theMzrObject,
                        std::string* fileName)
{
    // The arguments we look for/accept are
    // '-g NUMBER'       Sets the generation depth to NUMBER
    // '-f FILENAME'     Loads the inputfile FILENAME


    // Skip the command name.
    --argc; ++argv;

    bool filenameWasSeen ( false );

    //  Grab the arguments one by one.
    while (0 < argc)
    {
        std::string arg (*argv);
        argv++;
        argc--;

        if (arg == "-g" || arg == "--generation-depth")
        {
            std::string strDepth = utl::mustGetArg(argc, argv);
            int depth = utl::argMustBeNNInt (strDepth);
            theMzrObject.setGenerateDepth ( depth );
        }
        else if (arg == "-f" || arg == "--file")
        {
            filenameWasSeen = true;
            *fileName = utl::mustGetArg (argc, argv);
        }
        else if (arg == "-h" || arg == "--help")
        {
            displayHelpMessage();
            exit(0);
        }
        else throw utl::unkArgXcpt (arg);
    }

    if (!filenameWasSeen)
    {
        throw utl::badFileNameXcpt();
    }
}


void displayHelpMessage()
{
    using std::cout;
    using std::endl;
    
    cout << "This is a wrapper that loads a libmoleculizer object, attaches a file and quits." << std::endl;

}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                                                                          
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
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//              
//

#ifndef DEBUG_HH
#define DEBUG_HH

#include "utl/stdIncludes.hh"
#include "utl/xcpt.hh"

namespace mzr
{

template <typename DebugClassT>
class InteractiveModeManager
{
public:
typedef void (DebugClassT::*DebugClassFunctionPtr)(void);

InteractiveModeManager(DebugClassT* classInstanceToDebug)
:
ptrDebugClass( classInstanceToDebug )
{}

void
runInteractiveMode()
{

unsigned int result = 0;

while(true)
{
std::cout << "########################################\n";
std::cout << "## Debug Mode Options" << endl;
std::cout << "########################################\n";

std::cout << "0:\tQuit" << endl;
unsigned int funcNdx = 1;

for( typename std::vector<std::pair<std::string, DebugClassFunctionPtr> >::const_iterator iter = theDebugFunctions.begin();
iter != theDebugFunctions.end();
++iter, ++funcNdx)
{
std::cout<< funcNdx << ":\t" << iter->first << endl;
}


std::cout << "Input:\t";
std::cin >> result;
std::cout << std::endl;

if ( result == 0 )
{
break;
}

if (result > theDebugFunctions.size() )
{
std::cerr << "Error: Input '" << result << "' out of range." << endl;
std::cerr << "Acceptable input range:\t0-" << theDebugFunctions.size() << endl;
continue;
}
// Call the chosen functions.

try
{
(ptrDebugClass->*theDebugFunctions[result - 1].second)();
pause();
}
catch(utl::xcpt e)
{
std::cerr << "Exception caught in InteractiveModeManager." << endl;
e.warn();
std::cerr << "Continuing..." << endl;
}
}


}

bool
addFunction(std::string name, DebugClassFunctionPtr theFuncPtr)
{
theDebugFunctions.push_back( std::make_pair( name, theFuncPtr) );
return true;
}

private:
void
pause() const
{
// cout << "(Hit Enter)";
// std::cin.get();
return;
}

DebugClassT* ptrDebugClass;
std::vector<std::pair<std::string, DebugClassFunctionPtr> > theDebugFunctions;

};

}
#endif

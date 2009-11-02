//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 Nathan Addy
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


#ifndef AUX_HPP
#define AUX_HPP

#include <string>

enum RunMode { Full = 0,
               BoundedRun = 1,
               Iterations = 2 };

struct inputArgsStruct
{
    bool fileDefined;
    bool fileIsXml;
    bool fileIsRules;
    std::string fileName;

    bool hasOutputFile;
    std::string outputFile;

    int maxSpecies;
    int maxRxns;
    int numberIters;

    bool verbose;
    bool quiet;

    int runMode;

    inputArgsStruct()
        :
        fileDefined( false ),
        fileIsXml( false ),
        fileIsRules(false ),
        fileName( "" ),
        hasOutputFile( false ),
        outputFile( "" ),
        maxSpecies( -1 ),
        maxRxns( -1 ),
        verbose( false ),
        quiet( false ),
        runMode( Full )
    {}
};


int main(int argc, char* argv[]);
void printAllSpeciesStreams(mzr::moleculizer& refMolzer);
void printAllReactions( mzr::moleculizer& refMolzer);
void printStreamByName( mzr::moleculizer& refMolzer, const std::string& streamName);
void printStreamByTag( mzr::moleculizer& refMolzer, const std::string& streamName);
bool getUninitializedSpecies( const mzr::moleculizer& moleculizerRef, std::string& speciesName);
void processCommandLineArgs( int argc, char* argv[], inputArgsStruct& theInputArgs);
void printAllSpeciesByName(mzr::moleculizer& theMolzer);
void printAllSpeciesByID(mzr::moleculizer& theMolzer, std::string str);
void printAllPlexFamilies( mzr::moleculizer& theMolzer, std::string str);
void displayHelpAndExitProgram();

mzr::moleculizer::CachePosition
createBoundedNetwork(mzr::moleculizer& refMolzer, int maxSpec, int maxRxns);

void createFullNetwork(mzr::moleculizer& refMolzer);
void doNIterations(mzr::moleculizer& refMolzer, int number);
void runFullRunMoleculizer( mzr::moleculizer& mzr, const inputArgsStruct& inputArgsStruct);
void runIterationsMoleculizer( mzr::moleculizer& mzr, const inputArgsStruct& inputArgsStruct);
void runBoundedRunMoleculizer( mzr::moleculizer& mzr, const inputArgsStruct& inputArgsStruct);

#endif

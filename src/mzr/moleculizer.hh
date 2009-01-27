//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 The Molecular Sciences Institute.
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


#ifndef MOLECULIZER_H
#define MOLECULIZER_H

#include "utl/defs.hh"
#include "mzr/mzrException.hh"
#include "mzr/unit.hh"
#include "fnd/reactionNetworkDescription.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrReaction.hh"


namespace mzr
{
    class unitsMgr;
    
    // The main bulk of this class can be found in ReactionNetworkDescription.
    class moleculizer :
        public fnd::ReactionNetworkDescription<mzrSpecies, mzrReaction>
    {
        
    public:
        // Moleculizer only has basic constructors and destructors.  
        moleculizer( void );
        virtual ~moleculizer( void );
        
    public:
        
        void generateCompleteNetwork();
        void generateCompleteNetwork(long maxNumSpecies, long maxNumRxns = -1);
        
    public:

        void DEBUG_sayHello() const;
        
        void attachFileName( const std::string& aFileName );
        void attachString( const std::string& documentAsString );
        void attachDocument( xmlpp::Document* pDoc );

        bool getModelHasBeenLoaded() const;
        
    public:
        
        int getGenerationDepth( void ) const;
        void setGenerateDepth( unsigned int generateDepth );
        
        void setRateExtrapolation( bool rateExtrapolation );
        bool getRateExtrapolation() const;
        
        int 
        getNumberOfPlexFamilies() const;

        //////////////////////////////////////////////////
        // 
        // Functions for working with species streams
        //
        //////////////////////////////////////////////////

        void 
        getSpeciesStreams( std::vector<std::string>& speciesStreamNames) const;
        
        
        int 
        getNumberOfSpeciesInSpeciesStream(const std::string& streamName) const;
        
        void 
        getSpeciesInSpeciesStream(const std::string& streamName,
                                  std::vector<const mzr::mzrSpecies*>& speciesVector) const;

        void 
        getSpeciesInSpeciesStream(const std::string& streamName,
                                  std::vector<mzr::mzrSpecies*>& speciesVector);

        
        
    public:
        const mzrSpecies* 
        getSpeciesWithName( const std::string& speciesName ) throw( mzr::IllegalNameXcpt );
        
    public:
        xmlpp::Document*
        makeDomOutput( bool verbose ) throw( std::exception );

        void writeOutputFile( const std::string& fileName, bool verbose = false);
        
    protected:
        void setModelHasBeenLoaded( bool value );
        
        void insertGeneratedNetwork( xmlpp::Element* generatedNetworkElement, bool verbose );
        
        
    public:
        void recordUserNameToGeneratedNamePair( const std::string& userName,
                                                const std::string& genName )
        {
            userNameToSpeciesIDChart.insert( std::make_pair( userName, genName ) );
        }
        
        bool
        nameIsUserName(const std::string& possibleUserName) const
        {
            return (userNameToSpeciesIDChart.find( possibleUserName) != userNameToSpeciesIDChart.end());
        }
        
        std::string
        convertUserNameToGeneratedName(const std::string& possibleUserName) const 
            throw( utl::xcpt )
        {
            std::map<std::string, std::string>::const_iterator iter( userNameToSpeciesIDChart.find(possibleUserName) );
            if (iter == userNameToSpeciesIDChart.end() ) throw mzr::unknownUserNameXcpt( "Error, UserName doesn't exist and thus cannot be converted to a user name.");
            
            return iter->second;
        }
        
    protected:
        void
        constructorPrelude( void );
        
        void
        verifyInput( const xmlpp::Element* const pRootElt,
                     const xmlpp::Element* const pModelElt,
                     const xmlpp::Element* const pStreamElt ) const
            throw( std::exception );
        
        void
        constructorCore( xmlpp::Element* pRootElt,
                         xmlpp::Element* pModelElt,
                         xmlpp::Element* pStreamElt )
            throw( std::exception );
        
        void
        moleculizerParseDomInput( xmlpp::Element* pRootElt,
                                  xmlpp::Element* pModelElt,
                                  xmlpp::Element* pStreamElt )
            throw( std::exception );
        
        std::map<std::string, std::string> userNameToSpeciesIDChart;
        
    private:
        ////////////////////////////////////////////////
        // Units loaded by the user, waiting for destruction.
        //
        // This class is the manager for units, and the place that new units
        // can be installed.  It's public because units need to get to
        // each other.
        unitsMgr* pUserUnits;
        
        // Codes the input capabilities of moleculizer, including its parsing
        // routine.
        inputCapabilities inputCap;
        
        static int DEFAULT_GENERATION_DEPTH;
        
        
        bool modelLoaded;
        bool extrapolationEnabled;

        // Now we store a copy of the parser, so that people can get a copy of the rules, at any time.
        xmlpp::DomParser theParser;
    };
    
}

#endif

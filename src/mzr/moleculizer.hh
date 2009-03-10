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
#include "mzr/pythonRulesManager.hh"


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

        CachePosition
        generateCompleteNetwork(long maxNumSpecies, long maxNumRxns = -1);
        
    public:

        void loadXmlFileName( const std::string& aFileName );
        void loadXmlString( const std::string& documentAsString );

        void loadCommonRulesFileName(const std::string& aFileName );
        void loadCommonRulesString(const std::string& commonRulesAsString);

        void loadParsedDocument( xmlpp::Document* pDoc );



        bool getModelHasBeenLoaded() const;


    public:
      void addParameterStatement(const std::string& statement );
      void addModificationStatement( std::string& statement);
      void addMolsStatement( std::string& statement);
      void addAllostericPlexStatement( std::string& statement);
      void addAllostericOmniStatement( std::string& statement);
      void addDimerizationGenStatement( std::string& statement);
      void addOmniGenStatement( std::string& statement);
      void addUniMolGenStatement( std::string& statement);
      void addSpeciesStreamStatement( std::string& statement);
        
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

        int getNumberOfSpeciesStreams() const;

      bool speciesWithTagIsInSpeciesStream(const std::string speciesTag, const std::string& speciesStream ) const;
      bool speciesWithUniqueIDIsInSpeciesStream(const std::string speciesTag, const std::string& speciesStream ) const;
                
        int 
        getNumberOfSpeciesInSpeciesStream(const std::string& streamName) const;
        
        void 
        getSpeciesInSpeciesStream(const std::string& streamName,
                                  std::vector<const mzr::mzrSpecies*>& speciesVector) const;

    public:
        const mzrSpecies* 
        getSpeciesWithUniqueID( SpeciesIDCref uniqueID ) throw( mzr::IllegalNameXcpt );
        
    public:
        xmlpp::Document*
        makeDomOutput( bool verbose ) throw( std::exception );

        xmlpp::Document*
        makeDomOutput( bool verboseXML, CachePosition networkSizeRange ) throw( std::exception );

    void 
    loadGeneratedNetwork( xmlpp::Element* pGeneratedNetworkElmt);

        void writeOutputFile( const std::string& fileName, bool verbose = false);
        void writeOutputFile( const std::string& fileName, bool verbose, CachePosition pos);
        
    protected:
        void setModelHasBeenLoaded( bool value );
        
        void insertGeneratedNetwork( xmlpp::Element* generatedNetworkElt, CachePosition pos, bool verbose );
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
        convertUserNameToSpeciesID(const std::string& possibleUserName) const 
            throw( utl::xcpt )
        {
            std::map<std::string, std::string>::const_iterator iter( userNameToSpeciesIDChart.find(possibleUserName) );
            if (iter == userNameToSpeciesIDChart.end() ) throw mzr::unknownUserNameXcpt( "Error, UserName doesn't exist and thus cannot be converted to a user name.");
            
            return iter->second;
        }

        std::string
        convertUserNameToTaggedName( const std::string& possibleUserName) const
            throw( utl::xcpt)
        {

            std::map<std::string, std::string>::const_iterator iter( userNameToSpeciesIDChart.find(possibleUserName) );
            if (iter == userNameToSpeciesIDChart.end() ) throw mzr::unknownUserNameXcpt( "Error, UserName doesn't exist and thus cannot be converted to a user name.");
            
            std::string speciesID(iter->second);

            return convertSpeciesIDToSpeciesTag(speciesID);
        }

        void 
        getUserNames(std::vector<std::string>& refVector) const;

        int getNumberOfDefinedModifications() const;
        int getNumberOfDefinedMols() const;
        int getNumberOfDefinedRules() const;

        int getNumberOfDimerReactionRules() const;
        int getNumberOfOmniGenReactionRules() const;
        int getNumberOfUniMolReactionRules() const;
        int getNumberOfReactionRules() const;


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

      PythonRulesManager rulesManager;

        // Now we store a copy of the parser, so that people can get a copy of the rules, at any time.
        xmlpp::DomParser theParser;

    };

    class restoreGeneratedSpecies
    {
      
    public:
        restoreGeneratedSpecies( mzr::moleculizer& refMolzer )
            :
            rMolzer( refMolzer )
        {}
        
        void operator()( const xmlpp::Node* pNode )
        {


            const xmlpp::Element* pElement = dynamic_cast<const xmlpp::Element*>( pNode );

            assert(pElement);


            std::string uniqueID = utl::dom::mustGetAttrString( pElement, "unique-id" );
            std::string isExpanded = utl::dom::mustGetAttrString( pElement, "expanded");

            const mzrSpecies* pSpecies = rMolzer.getSpeciesWithUniqueID( uniqueID );

            std::string tag = pSpecies->getTag();
            
            if ( isExpanded == "true" )
            {
                rMolzer.incrementNetworkBySpeciesTag( tag );
            }
        }

    private:
        moleculizer& rMolzer;
        
    };
}

#endif

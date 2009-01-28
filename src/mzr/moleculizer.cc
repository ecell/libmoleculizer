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

#include <libxml++/libxml++.h>
#include <functional>

#include "utl/defs.hh"
#include "utl/utility.hh"
#include "utl/utlXcpt.hh"
#include "utl/dom.hh"
#include "utl/linearHash.hh"

#include "mzr/moleculizer.hh"
#include "mzr/mzrException.hh"

#include "mzr/mzrUnit.hh"

#include "mzr/unitsMgr.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/inputCapTest.hh"


namespace mzr
{
    // For getting each unit to do its part of parsing input.
    class unitParseDomInput :
        public std::unary_function<unit*, void>
    {
        xmlpp::Element* pRootElt;
        xmlpp::Element* pModelElt;
        xmlpp::Element* pStreamsElt;
        
    public:
        unitParseDomInput( xmlpp::Element* pRootElement,
                           xmlpp::Element* pModelElement,
                           xmlpp::Element* pStreamsElement ) :
            pRootElt( pRootElement ),
            pModelElt( pModelElement ),
            pStreamsElt( pStreamsElement )
        {}
        void
        operator()( unit* pUnit ) const throw( std::exception )
        {
            pUnit->parseDomInput( pRootElt,
                                  pModelElt,
                                  pStreamsElt );
        }
    };
    
    class prepareUnitToRun :
        public std::unary_function<unit*, void>
    {
        xmlpp::Element* pRootElt;
        xmlpp::Element* pModelElt;
        xmlpp::Element* pStreamsElt;
        
    public:
        prepareUnitToRun( xmlpp::Element* pRootElement,
                          xmlpp::Element* pModelElement,
                          xmlpp::Element* pStreamsElement ) :
            pRootElt( pRootElement ),
            pModelElt( pModelElement ),
            pStreamsElt( pStreamsElement )
        {
        }
        
        void
        operator()( unit* pUnit ) const
            throw( std::exception )
        {
            pUnit->prepareToRun( pRootElt,
                                 pModelElt,
                                 pStreamsElt );
        }
    };
    
    void
    moleculizer::constructorPrelude( void )
    {
        // Add up the input capabilities of the units.
        pUserUnits->unionInputCaps( inputCap );
        
        // Set the default generation depth.
        setGenerateDepth( moleculizer::DEFAULT_GENERATION_DEPTH );
    }
    
    void
    moleculizer::constructorCore( xmlpp::Element* pRootElement,
                                  xmlpp::Element* pModelElement,
                                  xmlpp::Element* pStreamsElement )
        throw( std::exception )
    {
        // Verify that the input can all be sucessfully handled by
        // various units.
        
        try
        {
            verifyInput( pRootElement,
                         pModelElement,
                         pStreamsElement );
        }
        catch ( utl::xcpt e )
        {
            e.warn();
        }
        
        moleculizerParseDomInput( pRootElement,
                                  pModelElement,
                                  pStreamsElement );
        
        std::for_each( pUserUnits->begin(),
                       pUserUnits->end(),
                       unitParseDomInput( pRootElement,
                                          pModelElement,
                                          pStreamsElement ) );
    }
    
    void
    moleculizer::moleculizerParseDomInput( xmlpp::Element* pRootElt,
                                           xmlpp::Element* pModelElt,
                                           xmlpp::Element* pStreamElt ) throw( std::exception )
    {
        ; // Do nothing for now.
    }
    
    moleculizer::moleculizer( void )
        :
        modelLoaded( false ),
        extrapolationEnabled( false ),
        theParser()
    {
        theParser.set_validate( false );

        pUserUnits = new unitsMgr( *this );
        
        // Now just does the "input capabilities" thing.
        constructorPrelude();
    }
    
    moleculizer::~moleculizer( void )
    {
        delete pUserUnits;
    }
    
    const mzrSpecies*
    moleculizer::getSpeciesWithName( const std::string& speciesName )
        throw( mzr::IllegalNameXcpt )
    {
        
        // First we check to see if it is a mangled name.
        mzrSpecies* theMzrSpecies;
        try
        {
            theMzrSpecies = findSpecies( speciesName );
            
            return theMzrSpecies;
        }
        catch ( fnd::NoSuchSpeciesXcpt e )
        {}
        
        // If not, we see if it is a user-name
        if (nameIsUserName( speciesName ) ) 
        {
            std::string mangledName = convertUserNameToGeneratedName( speciesName );
            return findSpecies( mangledName );
        }
        
        
        try
        {
            theMzrSpecies = pUserUnits->pNmrUnit->constructSpeciesFromName( speciesName );
            // Does mzrSpecies need to be expanded here?
            // theMzrSpecies->expandReactionNetwork();
            
            return theMzrSpecies;
        }
        catch ( nmr::IllegalNameXcpt xcpt )
        {
            throw mzr::IllegalNameXcpt( xcpt.getName() );
        }
    }
    
    void moleculizer::attachFileName( const std::string& filename )
    {
        theParser.parse_file( filename );
        if ( !theParser ) throw utl::dom::noDocumentParsedXcpt();
        this->attachDocument( theParser.get_document() );
    }
    
    void moleculizer::attachString( const std::string& documentAsString )
    {
        theParser.parse_memory( documentAsString );
        if ( !theParser ) throw utl::dom::noDocumentParsedXcpt();
        
        this->attachDocument( theParser.get_document() );
    }
    
    bool
    moleculizer::getModelHasBeenLoaded() const
    {
        return modelLoaded;
    }
    
    void
    moleculizer::setModelHasBeenLoaded( bool value )
    {
        modelLoaded = value;
    }

    void 
    moleculizer::DEBUG_sayHello() const
    {
        std::cout << "Hello from Libmoleculizer!!!" << std::endl;
    }

    void
    moleculizer::attachDocument( xmlpp::Document* pDoc )
    {
        if ( getModelHasBeenLoaded() )
        {
            throw utl::modelAlreadyLoadedXcpt();
        }
        else
        {
            setModelHasBeenLoaded( true );
        }
        
        // Get the basic framework of the document.
        xmlpp::Element* pRootElement
            = pDoc->get_root_node();
        
        // Get the unique model element.
        xmlpp::Element* pModelElement
            = utl::dom::mustGetUniqueChild( pRootElement,
                                            eltName::model );
        
        xmlpp::Element* pStreamsElement = utl::dom::getOptionalChild( pRootElement,
                                                                      eltName::streams );
        
        // Extract model info.
        constructorCore( pRootElement,
                         pModelElement,
                         pStreamsElement );
        
        // Have each unit do its prepareToRun thing.
        std::for_each( pUserUnits->begin(),
                       pUserUnits->end(),
                       prepareUnitToRun( pRootElement,
                                         pModelElement,
                                         pStreamsElement ) );
    }
    
    
    void
    moleculizer::setGenerateDepth( unsigned int generateDepth )
    {
        mzrUnit& rMzrUnit = * ( pUserUnits->pMzrUnit );
        rMzrUnit.setGenerateDepth( generateDepth );
    }
    
    int moleculizer::getGenerationDepth() const
    {
        return pUserUnits->pMzrUnit->getGenerateDepth();
    }
    
    class unitInsertStateElements :
        public std::unary_function<unit*, void>
    {
        xmlpp::Element* pModelElt;
    public:
        unitInsertStateElements( xmlpp::Element* pModelElement ) :
            pModelElt( pModelElement )
        {}
        
        void
        operator()( unit* pUnit ) const throw( std::exception )
        {
            pUnit->insertStateElts( pModelElt );
        }
    };
    
    xmlpp::Document*
    moleculizer::makeDomOutput( bool verbose ) throw( std::exception )
    {
        xmlpp::Document* pDoc = new xmlpp::Document();
        
        // Create the moleculizer-state node.
        xmlpp::Element* pRootElt
            = pDoc->create_root_node( eltName::moleculizerState );

        // Copy in the original model and streams stuff...
        xmlpp::Document* originalDoc = theParser.get_document();
        xmlpp::Element* originalRoot = originalDoc->get_root_node();

        xmlpp::Element* pInputModelElt = utl::dom::mustGetUniqueChild( originalRoot, eltName::model);
        xmlpp::Element* pInputStreamsElt = utl::dom::getOptionalChild( originalRoot, eltName::streams);

        pRootElt->import_node( pInputModelElt );
        if (pInputStreamsElt) pRootElt->import_node( pInputStreamsElt );
        
        xmlpp::Element* generatedNetworkElt = pRootElt->add_child( eltName::generatedNetwork );
        xmlpp::Element* unitStatesElement = pRootElt->add_child( eltName::unitsStates );
        
        // Insert each of the state elements
        this->insertGeneratedNetwork( generatedNetworkElt, verbose );
        
        // Run through the units, letting each make its complete state
        // contribution, for now.
        std::for_each( pUserUnits->begin(),
                       pUserUnits->end(),
                       unitInsertStateElements( unitStatesElement ) );
        
        return pDoc;
    }

    void moleculizer::writeOutputFile( const std::string& fileName, bool verbose)
    {
        xmlpp::Document* outputDocument = makeDomOutput(verbose);

        outputDocument->write_to_file_formatted( fileName );

        delete outputDocument;
    }


    void moleculizer::insertGeneratedNetwork( xmlpp::Element* generatedNetworkElt, bool verbose )
    {
        xmlpp::Element* genSpecElt = generatedNetworkElt->add_child("generated-species");
        xmlpp::Element* genRxnsElt = generatedNetworkElt->add_child("generated-reactions");


        for( SpeciesCatalogCIter specIter = this->getSpeciesCatalog().begin();
             specIter != this->getSpeciesCatalog().end();
             ++specIter)
        {
            xmlpp::Element* newSpecElt = genSpecElt->add_child( "species");
            newSpecElt->set_attribute("tag", *(specIter->first));
            newSpecElt->set_attribute("unique-id", convertSpeciesTagToSpeciesID( *specIter->first));

            if( specIter->second->hasNotified() ) 
            {
                newSpecElt->set_attribute("expanded", "true" );
            }
            else
            {
                newSpecElt->set_attribute("expanded", "false" );
            }
        }

       
        for( ReactionList::const_iterator rxnIter = this->getReactionList().begin();
             rxnIter != this->getReactionList().end();
             ++rxnIter)
        {
            const mzr::mzrReaction* pRxn( *rxnIter );

            xmlpp::Element* newRxnElt = genSpecElt->add_child( "reaction");
            xmlpp::Element* newSubstratesElt = newRxnElt->add_child("substrates");
            xmlpp::Element* newProductsElt = newRxnElt->add_child("products");
            xmlpp::Element* newRateElt = newRxnElt->add_child("rate");

            // newRxnElt->set_attribute( "representation", pRxn->getName() );

            // process each substrate

            for( mzr::mzrReaction::multMap::const_iterator reactantIter = pRxn->getReactants().begin();
                 reactantIter != pRxn->getReactants().end();
                 ++reactantIter )
            {
                
                const mzr::mzrReaction::multMap::value_type& vt( *reactantIter );

                xmlpp::Element* newSubElt = newSubstratesElt->add_child("substrate");
                newSubElt->set_attribute("multiplicity", utl::stringify(vt.second) );
                newSubElt->set_attribute("tag", vt.first->getTag() );
                if(verbose)
                {
                    newSubElt->set_attribute("unique-id", vt.first->getName() );
                }
            }

            // Process each reactant
            for( mzr::mzrReaction::multMap::const_iterator productIter = pRxn->getProducts().begin();
                 productIter != pRxn->getProducts().end();
                 ++productIter )
            {
                const mzr::mzrReaction::multMap::value_type& vt( *productIter );

                xmlpp::Element* newProdElt = newProductsElt->add_child("product");
                newProdElt->set_attribute("multiplicity", utl::stringify(vt.second) );
                newProdElt->set_attribute("tag", vt.first->getTag() );
                if(verbose)
                {
                    newProdElt->set_attribute("unique-id", vt.first->getName() );
                }
            }

            newRateElt->set_attribute("value", utl::stringify(pRxn->getRate() ) );
        }
        
    }
    
    void
    moleculizer::verifyInput( xmlpp::Element const * const pRootElement,
                              xmlpp::Element const * const pModelElement,
                              xmlpp::Element const * const pStreamsElement ) const
        throw( std::exception )
    {
        // Verify that every element in the model section is handled by some
        // unit or another.
        //
        // I discover very late in the game that this code incorrectly tries
        // to convert comments to elements.
        xmlpp::Node::NodeList modelContentNodes
            = pModelElement->get_children();

        xmlpp::Node::NodeList::iterator iUnhandledModelContent
            = std::find_if( modelContentNodes.begin(),
                            modelContentNodes.end(),
                            modelNodeNotInCap( inputCap ) );

        if ( modelContentNodes.end() != iUnhandledModelContent )
            throw unhandledModelContentXcpt( *iUnhandledModelContent );
        
        // Get the reaction-gens node, which contains kinds of reaction generators
        // introduced by units.
        xmlpp::Element* pReactionGensElt
            = utl::dom::mustGetUniqueChild( pModelElement,
                                            eltName::reactionGens );
        
        // Verify that every kind of reaction generator is handled by some
        // unit or another.
        xmlpp::Node::NodeList reactionGensContentNodes
            = pReactionGensElt->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledReactionGenContent
            = std::find_if( reactionGensContentNodes.begin(),
                            reactionGensContentNodes.end(),
                            reactionGenNotInCap( inputCap ) );
        
        if ( reactionGensContentNodes.end() !=
             iUnhandledReactionGenContent )
            throw unhandledReactionGenXcpt( *iUnhandledReactionGenContent );
        
        // Get the explicit species model node, which can contain kinds of explicit
        // species introduced by units.
        xmlpp::Element* pExplicitSpeciesElt
            = utl::dom::mustGetUniqueChild( pModelElement,
                                            eltName::explicitSpecies );
        
        // Verify that every kind of explicit species is handled by some
        // unit or another.
        xmlpp::Node::NodeList explicitSpeciesContentNodes
            = pExplicitSpeciesElt->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledExplicitSpeciesContent
            = std::find_if( explicitSpeciesContentNodes.begin(),
                            explicitSpeciesContentNodes.end(),
                            explicitSpeciesNodeNotInCap( inputCap ) );
        
        if ( explicitSpeciesContentNodes.end() != iUnhandledExplicitSpeciesContent )
            throw unhandledExplicitSpeciesContentXcpt( *iUnhandledExplicitSpeciesContent );
        
    }

    void 
    moleculizer::generateCompleteNetwork() 
    {
        if ( ! this->getModelHasBeenLoaded() ) throw ModelNotLoadedXcpt("moleculizer::generateCompleteNetwork");

        // This function will generate the entire network.  To be used primarily 
        bool nodesUnexpanded = true;
        
        while( nodesUnexpanded)
        {
            for( SpeciesCatalogCIter speciesIter = this->getSpeciesCatalog().begin();
                 speciesIter != this->getSpeciesCatalog().end();
                 ++speciesIter)
            {
                speciesIter->second->expandReactionNetwork();
            }
            
            nodesUnexpanded = false;

            for( SpeciesCatalogCIter speciesIter = this->getSpeciesCatalog().begin();
                 speciesIter != this->getSpeciesCatalog().end();
                 ++speciesIter)
            {
                if (speciesIter->second->hasNotified() )
                {
                    nodesUnexpanded = true;
                    break;
                }
            }
        }
    }

    void moleculizer::generateCompleteNetwork(long maxNumSpecies, long maxNumRxns)
    {
        if ( ! this->getModelHasBeenLoaded() ) throw ModelNotLoadedXcpt("moleculizer::generateCompleteNetwork");

        if (maxNumSpecies < 0) maxNumSpecies = std::numeric_limits<long>::max();
        if (maxNumRxns < 0) maxNumRxns = std::numeric_limits<long>::max();


        // This loop has some goofy logic.  Is there a better way to do this?  Seems like there 
        // has to be...
        
        const std::string NULLSTRING("");
        while( maxNumSpecies > getTotalNumberSpecies() && 
               maxNumRxns > getTotalNumberReactions() )
        {
            std::string speciesToInc("");


            for( SpeciesCatalogCIter speciesIter = this->getSpeciesCatalog().begin();
                 speciesIter != this->getSpeciesCatalog().end();
                 ++speciesIter)
            {
                if (! speciesIter->second->hasNotified() )
                {
                    speciesToInc = *(speciesIter->first);
                    break;  // Out of the for each loop
                }
            }

            if (speciesToInc == NULLSTRING )
            {
                // From this we may imply that all species have been notified and hence we return.
                return;
            }
            else
            {
                this->incrementNetworkBySpeciesTag( speciesToInc );
            }
            
        }
    }
    
    bool moleculizer::getRateExtrapolation( void ) const
    {
        return extrapolationEnabled;
    }
    
    void moleculizer::setRateExtrapolation( bool rateExtrapolation ) 
    {

        if ( getModelHasBeenLoaded() )
        {
            throw utl::modelAlreadyLoadedXcpt();
        }

        extrapolationEnabled = rateExtrapolation;
    }
    
    int 
    moleculizer::getNumberOfSpeciesInSpeciesStream(const std::string& streamName) const 
    {
        // TODO
        return 0;
    }

    void 
    moleculizer::getSpeciesInSpeciesStream(const std::string& streamName,
                                           std::vector<const mzr::mzrSpecies*>& speciesVector) const
    {
        fnd::dumpable<fnd::basicDumpable::dumpArg>* ptrDumpable = 
            pUserUnits->pMzrUnit->mustFindDumpable( streamName );
        
        
        mzr::multiSpeciesDumpable<plx::mzrPlexSpecies>* streamPtr;
        
        streamPtr = dynamic_cast< mzr::multiSpeciesDumpable<plx::mzrPlexSpecies>* >( ptrDumpable);
        
        const std::vector<const plx::mzrPlexSpecies*>* theSpecies = streamPtr->getSpeciesInMultiSpeciesStream();
        
        
        if (!streamPtr)
        {
            std::cerr << "Error finding dumpable " << streamName << std::endl;
            exit(1);
        }
        
        for( std::vector<const plx::mzrPlexSpecies*>::const_iterator ssIter = theSpecies->begin();
             ssIter != theSpecies->end();
             ++ssIter)
        {
            speciesVector.push_back( *ssIter );
        }
    }

    void 
    moleculizer::getSpeciesInSpeciesStream(const std::string& streamName,
                                           std::vector<mzr::mzrSpecies*>& speciesVector)
    {
        // Do nothing for now.
    }
    
    void 
    moleculizer::getSpeciesStreams( std::vector<std::string>& speciesStreamNames) const
    {
        pUserUnits->pMzrUnit->getListOfDumpables( speciesStreamNames );
    }
    
    int moleculizer::getNumberOfPlexFamilies() const
    {
        return pUserUnits->pPlexUnit->familyCount();
    }
    
    
    
    int moleculizer::DEFAULT_GENERATION_DEPTH = 1;
    
}

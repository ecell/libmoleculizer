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
#include <iterator>

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
#include "dimer/dimerUnit.hh"
#include "ftr/ftrUnit.hh"


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
    moleculizer::getSpeciesWithUniqueID( SpeciesIDCref uniqueid )
        throw( mzr::IllegalNameXcpt )
    {
        
        // First we check to see if it is a mangled name.
        mzrSpecies* theMzrSpecies;
        try
        {
            std::string tag = convertSpeciesIDToSpeciesTag( uniqueid );
            theMzrSpecies = findSpecies( tag );
            return theMzrSpecies;
        }
        catch ( fnd::NoSuchSpeciesXcpt e )
        {}
        
        try
        {
            theMzrSpecies = pUserUnits->pNmrUnit->constructSpeciesFromName( uniqueid );
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

        // Check to see if it has a "generated-network" tag, which would mean
        // it is an output mzr file being read in.
        xmlpp::Element* pGeneratedNetworkElmt = utl::dom::getOptionalChild( pRootElement,
                                                                            "generated-network");
        if (pGeneratedNetworkElmt)
        {
            this->loadGeneratedNetwork( pGeneratedNetworkElmt );
        }
    }

    void 
    moleculizer::loadGeneratedNetwork( xmlpp::Element* pGeneratedNetworkElmt)
    {

        xmlpp::Element* pGeneratedSpeciesElmt = \
            utl::dom::mustGetUniqueChild( pGeneratedNetworkElmt, "generated-species" ); 

        xmlpp::Node::NodeList genSpecNodes = pGeneratedSpeciesElmt->get_children("species");

        std::for_each( genSpecNodes.begin(),
                       genSpecNodes.end(),
                       restoreGeneratedSpecies( *this ) );

        return;

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
    moleculizer::makeDomOutput( bool verboseXML ) throw( std::exception )
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
        this->insertGeneratedNetwork( generatedNetworkElt, verboseXML );
        
        // Run through the units, letting each make its complete state
        // contribution, for now.
        std::for_each( pUserUnits->begin(),
                       pUserUnits->end(),
                       unitInsertStateElements( unitStatesElement ) );
        
        return pDoc;
    }

    xmlpp::Document*
    moleculizer::makeDomOutput( bool verboseXML, CachePosition pos ) throw( std::exception )
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
        this->insertGeneratedNetwork( generatedNetworkElt, pos, verboseXML );
        
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

    void moleculizer::writeOutputFile( const std::string& fileName, bool verbose, CachePosition pos)
    {
        xmlpp::Document* outputDocument = makeDomOutput(verbose, pos);

        outputDocument->write_to_file_formatted( fileName );

        delete outputDocument;
    }


    void moleculizer::insertGeneratedNetwork( xmlpp::Element* generatedNetworkElt, bool verbose )
    {
        xmlpp::Element* genSpecElt = generatedNetworkElt->add_child("generated-species");
        xmlpp::Element* genRxnsElt = generatedNetworkElt->add_child("generated-reactions");


        // Insert each of the generated species into the network.
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

        // Insert each of the generated reactions into the network.
        for( ReactionList::const_iterator rxnIter = this->getReactionList().begin();
             rxnIter != this->getReactionList().end();
             ++rxnIter)
        {
            const mzr::mzrReaction* pRxn( *rxnIter );

            xmlpp::Element* newRxnElt = genRxnsElt->add_child( "reaction");
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


    void moleculizer::insertGeneratedNetwork( xmlpp::Element* generatedNetworkElt, CachePosition pos, bool verbose )
    {
        xmlpp::Element* genSpecElt = generatedNetworkElt->add_child("generated-species");
        xmlpp::Element* genRxnsElt = generatedNetworkElt->add_child("generated-reactions");


        // Insert each of the generated species into the network.

        for( SpeciesListCIter specIter = this->getDeltaSpeciesList().begin();
             specIter != pos.first;
             ++specIter )
        {
            xmlpp::Element* newSpecElt = genSpecElt->add_child( "species");
            
            newSpecElt->set_attribute("tag", (*specIter)->getTag());
            newSpecElt->set_attribute("unique-id", convertSpeciesTagToSpeciesID( (*specIter)->getTag()));

            if( (*specIter)->hasNotified() ) 
            {
                newSpecElt->set_attribute("expanded", "true" );
            }
            else
            {
                newSpecElt->set_attribute("expanded", "false" );
            }
        }

        // Insert each of the generated reactions into the network.
        for( ReactionList::const_iterator rxnIter = this->getDeltaReactionList().begin();
             rxnIter != pos.second;
             ++rxnIter)
        {
            const mzr::mzrReaction* pRxn( *rxnIter );

            xmlpp::Element* newRxnElt = genRxnsElt->add_child( "reaction");
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

        // This function will generate the entire network.  
        
        // This is a rewrite, because the previous function seemed like it might not be working.


        bool foundUnexpandedSpecies = false;
        mzrSpecies* ptrUnexpandedSpecies = NULL;

        while(true)
        {
            foundUnexpandedSpecies = false;

            for( SpeciesCatalogCIter speciesIter = this->getSpeciesCatalog().begin();
                 speciesIter != this->getSpeciesCatalog().end();
                 ++speciesIter)
            {
                if ( !speciesIter->second->hasNotified() )
                {
                    foundUnexpandedSpecies = true;
                    ptrUnexpandedSpecies = speciesIter->second;
                    break;
                }
            }

            if (foundUnexpandedSpecies)
            {
                // Expand and continue...
                ptrUnexpandedSpecies->expandReactionNetwork();
            }
            else
            {
                // GeneratedCompleteNetwork
                break;
            }
        }

    }

    
    moleculizer::CachePosition
    moleculizer::generateCompleteNetwork(long maxNumSpecies, long maxNumRxns)
    {
        // Note when reading this function:
        // Lists have the property that no iterators are ever invalidated.  It turns out (GOTCHA!)
        // this means that when you define an iterator as equal to theList.end(), this iterator will
        // always point to the end(), no matter what is added to the list in the interim.

        // Therefore, this function works by getting the end() and then producing end()--.  This iterator 
        // will always point to "the last element in the list".  Then I ++ it at various places to get back what
        // would have been the end() previously.  Again, the decrementing and incrementing is so that I can get a 
        // never invalidated iterator that has the same value as end() but an absolute position.  
        // 
        //It is probably possible and desirable to replace this with end() and then increment it.  


        // I don't know how to enforce this at the moment, but it is important that 
        // neither delta cache has been cleared prior to this function running.  

        if ( ! this->getModelHasBeenLoaded() ) throw ModelNotLoadedXcpt("moleculizer::generateCompleteNetwork");

        if (maxNumSpecies < 0) maxNumSpecies = std::numeric_limits<long>::max();
        if (maxNumRxns < 0) maxNumRxns = std::numeric_limits<long>::max();

        // This isn't really necessary, but it should be true anyways (there is nothing to expand otherwise), and gets
        // us away from the concern regarding lists explicated at the beginning of this function.
        assert( this->getTotalNumberSpecies() > 0 && this->getTotalNumberReactions() > 0);

        if ( (long) this->getTotalNumberSpecies() >= maxNumSpecies || (long) this->getTotalNumberReactions() >= maxNumRxns )
        {
            // We are already too big to know what to do with ourselves.  Throw an exception.
            throw utl::xcpt("Error in moleculizer::generateCompleteNetwork(long maxNumSpecies, long maxNumRxns).  Network is a priori too large for parameters.");
        }

        SpeciesListIter specCacheMaxIter = theDeltaSpeciesList.end();
        ReactionListIter rxnCacheMaxIter = theDeltaReactionList.end();
        
        specCacheMaxIter--;
        rxnCacheMaxIter--;

        bool expandedOne = false;

        // We start out good -- the reaction network is not too big, because we passed the precondition
        while( true )
        {
            expandedOne = false;

            for( SpeciesCatalogCIter speciesIter = this->getSpeciesCatalog().begin();
                 speciesIter != this->getSpeciesCatalog().end();
                 ++speciesIter)
            {
                if ( !speciesIter->second->hasNotified() )
                {
                    expandedOne = true;
                    speciesIter->second->expandReactionNetwork();
                    break;
                }
            }

            // The network is too big, since we came into the loop good, return that value of cached stuff.
            if ( (long) getTotalNumberSpecies() > maxNumSpecies || (long) getTotalNumberReactions() > maxNumRxns )
            {
                return std::make_pair( ++specCacheMaxIter, ++rxnCacheMaxIter);
            }

            specCacheMaxIter = theDeltaSpeciesList.end();
            specCacheMaxIter--;
            rxnCacheMaxIter = theDeltaReactionList.end();
            rxnCacheMaxIter--;

            if (!expandedOne)
            {
                return std::make_pair( theDeltaSpeciesList.end(), theDeltaReactionList.end() );
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

    int moleculizer::getNumberOfSpeciesStreams() const
    {
        return pUserUnits->pMzrUnit->getNumberOfDumpables();
    }
    
    int moleculizer::getNumberOfPlexFamilies() const
    {
        return pUserUnits->pPlexUnit->familyCount();
    }


    int moleculizer::getNumberOfDefinedModifications() const
    {
        return pUserUnits->pMolUnit->getNumberKnownMods();
    }

    int moleculizer::getNumberOfDefinedMols() const
    {
        return pUserUnits->pMolUnit->getNumberKnownMols();
    }

    int moleculizer::getNumberOfDefinedRules() const
    {
        return getNumberOfDimerReactionRules() + 
            getNumberOfOmniGenReactionRules() + 
            getNumberOfUniMolReactionRules();
    }

    int moleculizer::getNumberOfDimerReactionRules() const
    {
        return pUserUnits->pDimerUnit->getNumberOfDimerDecompRules();
    }

    int moleculizer::getNumberOfOmniGenReactionRules() const
    {
        return pUserUnits->pFtrUnit->getNumberOfOmniGenRules();
    }
    int moleculizer::getNumberOfUniMolReactionRules() const
    {
        return pUserUnits->pFtrUnit->getNumberOfUniMolGenRules();
    }
    
    int moleculizer::DEFAULT_GENERATION_DEPTH = 1;
    
}

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

#include <libxml++/libxml++.h>
#include <boost/foreach.hpp>
#include <functional>

#include "utl/defs.hh"
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
        unitParseDomInput (xmlpp::Element* pRootElement,
                           xmlpp::Element* pModelElement,
                           xmlpp::Element* pStreamsElement) :
                pRootElt (pRootElement),
                pModelElt (pModelElement),
                pStreamsElt ( pStreamsElement)
        {}
        void
        operator() (unit* pUnit) const throw (std::exception)
        {
            pUnit->parseDomInput (pRootElt,
                                  pModelElt,
                                  pStreamsElt);
        }
    };

    class prepareUnitToRun :
                public std::unary_function<unit*, void>
    {
        xmlpp::Element* pRootElt;
        xmlpp::Element* pModelElt;
        xmlpp::Element* pStreamsElt;

    public:
        prepareUnitToRun (xmlpp::Element* pRootElement,
                          xmlpp::Element* pModelElement,
                          xmlpp::Element* pStreamsElement) :
                pRootElt (pRootElement),
                pModelElt (pModelElement),
                pStreamsElt (pStreamsElement)
        {
        }

        void
        operator() (unit* pUnit) const
        throw (std::exception)
        {
            pUnit->prepareToRun (pRootElt,
                                 pModelElt,
                                 pStreamsElt);
        }
    };

    void
    moleculizer::constructorPrelude (void)
    {
        // Add up the input capabilities of the units.
        pUserUnits->unionInputCaps (inputCap);

        // Set the default generation depth.
        setGenerateDepth( moleculizer::DEFAULT_GENERATION_DEPTH );
    }

    void
    moleculizer::constructorCore (xmlpp::Element* pRootElement,
                                  xmlpp::Element* pModelElement,
                                  xmlpp::Element* pStreamsElement)
    throw (std::exception)
    {
        // Verify that the input can all be sucessfully handled by
        // various units.

        try
        {
            verifyInput (pRootElement,
                         pModelElement,
                         pStreamsElement);
        }
        catch (utl::xcpt e)
        {
            e.warn();
        }

        std::for_each (pUserUnits->begin(),
                       pUserUnits->end(),
                       unitParseDomInput (pRootElement,
                                          pModelElement,
                                          pStreamsElement) );
    }

    moleculizer::moleculizer(void)
        :
        modelLoaded( false )
    {
        
        pUserUnits = new unitsMgr(*this);

        // Now just does the "input capabilities" thing.
        constructorPrelude();
    }

    moleculizer::~moleculizer (void)
    {
        delete pUserUnits;
    }

    mzrSpecies*
    moleculizer::getSpeciesWithName (const std::string& speciesName)
    throw ( mzr::IllegalNameXcpt )
    {

        mzrSpecies* theMzrSpecies;
        try
        {
            theMzrSpecies = findSpecies (speciesName);
//            theMzrSpecies->expandReactionNetwork();

            return theMzrSpecies;
        }
        catch (fnd::NoSuchSpeciesXcpt e)
        {}

        try
        {
            theMzrSpecies = pUserUnits->pNmrUnit->constructSpeciesFromName (speciesName);
            // Does mzrSpecies need to be expanded here?
            // theMzrSpecies->expandReactionNetwork();

            return theMzrSpecies;
        }
        catch (nmr::IllegalNameXcpt xcpt)
        {
            throw mzr::IllegalNameXcpt ( xcpt.getName() );
        }
    }

    std::string
    moleculizer::getRandomSpeciesName() const
    {
        // I use pointers here, to save space.
        // I should make this much more efficient.

        std::vector<mzrSpecies*> allSpecies;
        BOOST_FOREACH ( SpeciesCatalog::value_type vt, theSpeciesListCatalog)
        {
            allSpecies.push_back ( vt.second );
        }

        int random_index = rand() % allSpecies.size();

        return allSpecies[random_index]->getName();
    }


    void moleculizer::attachFileName (const std::string& filename)
    {
        xmlpp::DomParser parser;
        parser.set_validate (false);
        parser.parse_file ( filename );
        if (!parser) throw utl::dom::noDocumentParsedXcpt();

        this->attachDocument ( parser.get_document() );

    }

    void moleculizer::attachString ( const std::string& documentAsString )
    {
        xmlpp::DomParser parser;
        parser.set_validate ( false );
        parser.parse_memory ( documentAsString );
        if (!parser) throw utl::dom::noDocumentParsedXcpt();

        this->attachDocument ( parser.get_document() );
    }

    bool
    moleculizer::getModelHasBeenLoaded() const
    {
        return modelLoaded;
    }

    void
    moleculizer::setModelHasBeenLoaded (bool value)
    {
        modelLoaded = value;
    }

    void
    moleculizer::attachDocument (xmlpp::Document* pDoc)
    {
        // Note that this new pattern here will allow old style moleculizer
        // files, ie those that contain events, to be parsed and run.

        if (getModelHasBeenLoaded() )
        {
            throw utl::modelAlreadyLoadedXcpt();
        }
        else
        {
            setModelHasBeenLoaded (true);
        }

        // Get the basic framework of the document.
        xmlpp::Element* pRootElement
        = pDoc->get_root_node();

        // Get the unique model element.
        xmlpp::Element* pModelElement
        = utl::dom::mustGetUniqueChild (pRootElement,
                                        eltName::model);

        xmlpp::Element* pStreamsElement
        = utl::dom::getOptionalChild (pRootElement,
                                      eltName::streams);

        // Extract model info.
        constructorCore (pRootElement,
                         pModelElement,
                         pStreamsElement);

        // Have each unit do its prepareToRun thing.
        std::for_each (pUserUnits->begin(),
                       pUserUnits->end(),
                       prepareUnitToRun (pRootElement,
                                         pModelElement,
                                         pStreamsElement) );

    }


    void
    moleculizer::setGenerateDepth (unsigned int generateDepth)
    {
        mzrUnit& rMzrUnit = * (pUserUnits->pMzrUnit);
        rMzrUnit.setGenerateDepth (generateDepth);
    }

    class unitInsertStateElements :
                public std::unary_function<unit*, void>
    {
        xmlpp::Element* pModelElt;
    public:
        unitInsertStateElements (xmlpp::Element* pModelElement) :
                pModelElt (pModelElement)
        {}

        void
        operator() (unit* pUnit) const throw (std::exception)
        {
            pUnit->insertStateElts (pModelElt);
        }
    };

    // State dump, invoked in a dump-state event.
    xmlpp::Document*
    moleculizer::makeDomOutput (void) throw (std::exception)
    {
        xmlpp::Document* pDoc = new xmlpp::Document();

        // Create the moleculizer-state node.
        xmlpp::Element* pRootElt
        = pDoc->create_root_node (eltName::moleculizerState);

        // Add some mandatory elements.
        xmlpp::Element* pModelElt
        = pRootElt->add_child (eltName::model);
        pModelElt->add_child (eltName::unitsStates);
        pModelElt->add_child (eltName::explicitSpeciesTags);
        pModelElt->add_child (eltName::taggedSpecies);
        pModelElt->add_child (eltName::tagReactions);

        // The only reason for inserting these elements here seems
        // to have been (thought to be?) that various modules would insert
        // elements under them.
        pModelElt->add_child (eltName::time);
        pRootElt->add_child (eltName::streams);

        // Run through the units, letting each make its complete state
        // contribution, for now.
        std::for_each (pUserUnits->begin(),
                       pUserUnits->end(),
                       unitInsertStateElements (pRootElt) );

        return pDoc;
    }

    void
    moleculizer::verifyInput (xmlpp::Element const * const pRootElement,
                              xmlpp::Element const * const pModelElement,
                              xmlpp::Element const * const pStreamsElement) const
    throw (std::exception)
    {
        // Verify that every element in the model section is handled by some
        // unit or another.
        //
        // I discover very late in the game that this code incorrectly tries
        // to convert comments to elements.
        xmlpp::Node::NodeList modelContentNodes
        = pModelElement->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledModelContent
        = std::find_if (modelContentNodes.begin(),
                        modelContentNodes.end(),
                        modelNodeNotInCap (inputCap) );
        if (modelContentNodes.end() != iUnhandledModelContent)
            throw unhandledModelContentXcpt (*iUnhandledModelContent);

        // Get the reaction-gens node, which contains kinds of reaction generators
        // introduced by units.
        xmlpp::Element* pReactionGensElt
        = utl::dom::mustGetUniqueChild (pModelElement,
                                        eltName::reactionGens);

        // Verify that every kind of reaction generator is handled by some
        // unit or another.
        xmlpp::Node::NodeList reactionGensContentNodes
        = pReactionGensElt->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledReactionGenContent
        = std::find_if (reactionGensContentNodes.begin(),
                        reactionGensContentNodes.end(),
                        reactionGenNotInCap (inputCap) );
        if (reactionGensContentNodes.end() !=
                iUnhandledReactionGenContent)
            throw unhandledReactionGenXcpt (*iUnhandledReactionGenContent);

        // Get the explicit species model node, which can contain kinds of explicit
        // species introduced by units.
        xmlpp::Element* pExplicitSpeciesElt
        = utl::dom::mustGetUniqueChild (pModelElement,
                                        eltName::explicitSpecies);

        // Verify that every kind of explicit species is handled by some
        // unit or another.
        xmlpp::Node::NodeList explicitSpeciesContentNodes
        = pExplicitSpeciesElt->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledExplicitSpeciesContent
        = std::find_if (explicitSpeciesContentNodes.begin(),
                        explicitSpeciesContentNodes.end(),
                        explicitSpeciesNodeNotInCap (inputCap) );

        if (explicitSpeciesContentNodes.end() != iUnhandledExplicitSpeciesContent)
            throw unhandledExplicitSpeciesContentXcpt (*iUnhandledExplicitSpeciesContent);

    }

    void 
    moleculizer::recordPlexParameter( const std::string& plexName,
                                      const std::string& paramName,
                                      const double& paramValue)
    {
        if (paramName == "kD")
        {
            k_DChart.insert( std::make_pair( plexName, paramValue) );
        }
        else if (paramName == "radius")
        {
            radiusChart.insert( std::make_pair( plexName, paramValue ) );
        }
         
        else throw 666;
    }

    void 
    moleculizer::enableSpatialExtrapolation()
    {
        boost::function<void (const mzrReaction*)> cb1, cb2, cb3, cb4;

        cb1 = std::bind1st(std::mem_fun(&moleculizer::installRadiusForNewSpecies), this);
        cb2 = std::bind1st(std::mem_fun(&moleculizer::installKdForNewSpecies), this);
        cb3 = std::bind1st(std::mem_fun(&moleculizer::installKaForNewBinaryReaction), this);
        cb4 = std::bind1st(std::mem_fun(&moleculizer::installKForNewUnaryReaction), this);

        // Note that the order of these may matter. 
        addCallback( cb1 );
        addCallback( cb2 );
        addCallback( cb3 );
        addCallback( cb4 );
    }
    

    void 
    moleculizer::enableNonspatialExtrapolation() 
    {
        boost::function<void (const mzrReaction*)> cb1, cb2;

        cb1 = std::bind1st(std::mem_fun(&moleculizer::installKForNewBinaryReaction), this);
        cb2 = std::bind1st(std::mem_fun(&moleculizer::installKForNewUnaryReaction), this);

        // Note that the order of these may matter. 
        addCallback( cb1 );
        addCallback( cb2 );
    }

    int moleculizer::DEFAULT_GENERATION_DEPTH = 1;

}

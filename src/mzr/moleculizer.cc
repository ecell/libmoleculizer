/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008  Walter Lawrence (Larry) Lok.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#include <libxml++/libxml++.h>
#include "utl/noDocumentParsedXcpt.hh"
#include "utl/modelPreviouslyLoadedXcpt.hh"
#include "mzr/mzrException.hh"
#include <set>
#include <algorithm>
#include <fstream>
#include <ctime>
#include "utl/platform.hh"
#include "utl/linearHash.hh"

// Maybe this can be taken out too....
#include "utl/dom.hh"

#include "mzr/mzrSpecies.hh"
#include "mzr/mzrReaction.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/moleculizer.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/inputCapTest.hh"
#include "mzr/inputCapXcpt.hh"


namespace mzr
{
    // For getting each unit to do its part of parsing input.
    class unitParseDomInput :
        public std::unary_function<unit*, void>
    {
        xmlpp::Element* pRootElt;
        xmlpp::Element* pModelElt;

    public:
        unitParseDomInput(xmlpp::Element* pRootElement,
                          xmlpp::Element* pModelElement) :
            pRootElt(pRootElement),
            pModelElt(pModelElement)
        {}
        void
        operator()(unit* pUnit) const throw(std::exception)
        {
            pUnit->parseDomInput(pRootElt,
                                 pModelElt);
        }
    };

    class prepareUnitToRun :
        public std::unary_function<unit*, void>
    {
        xmlpp::Element* pRootElt;
        xmlpp::Element* pModelElt;

    public:
        prepareUnitToRun(xmlpp::Element* pRootElement,
                         xmlpp::Element* pModelElement) :
            pRootElt(pRootElement),
            pModelElt(pModelElement)
        {
        }

        void
        operator()(unit* pUnit) const
            throw(std::exception)
        {
            pUnit->prepareToRun(pRootElt,
                                pModelElt);
        }
    };

    void
    moleculizer::constructorPrelude(void)
    {
        // Add up the input capabilities of the units.
        pUserUnits->unionInputCaps(inputCap);
    }

    void
    moleculizer::constructorCore(xmlpp::Element* pRootElement,
                                 xmlpp::Element* pModelElement)
        throw(std::exception)
    {
        // Verify that the input can all be sucessfully handled by 
        // various units.

        try
        {
            verifyInput(pRootElement,
                        pModelElement);
        }
        catch(utl::xcpt e)
        {
            e.warn();
        }

        std::for_each(pUserUnits->begin(),
                      pUserUnits->end(),
                      unitParseDomInput(pRootElement,
                                        pModelElement));
    }

    moleculizer::moleculizer(void)         
        :
        modelLoaded( false )
    {
        
        pUserUnits = new unitsMgr(*this);

        // Now just does the "input capabilities" thing.
        constructorPrelude();
    }

    moleculizer::~moleculizer(void)
    {
        delete pUserUnits;
    }

  mzrSpecies*
  moleculizer::getSpeciesWithName(const std::string& speciesName) 
    throw( mzr::IllegalNameXcpt )
  {

    mzrSpecies* theMzrSpecies;
    try
      {
	theMzrSpecies = findSpecies(speciesName);

        // Don't know is this is strictly necessary...
        theMzrSpecies->expandReactionNetwork(); 

	return theMzrSpecies;
      }
    catch(fnd::NoSuchSpeciesXcpt xcpt)
      {
      }


    try
      {
        theMzrSpecies = pUserUnits->pNmrUnit->constructSpeciesFromName(speciesName);
        
        // Does mzrSpecies need to be expanded here?
        theMzrSpecies->expandReactionNetwork();
        
        return theMzrSpecies;
      }
    catch(nmr::IllegalNameXcpt xcpt)
      {
        throw mzr::IllegalNameXcpt( xcpt.getName() );
      }
  }


    void moleculizer::attachFileName(const std::string& filename)
    {
        xmlpp::DomParser parser;
        parser.set_validate(false);
        parser.parse_file( filename );
        if (!parser) throw utl::dom::noDocumentParsedXcpt();
            
        this->attachDocument( parser.get_document() );
    }

    void moleculizer::attachString( const std::string& documentAsString )
    {
        xmlpp::DomParser parser;
        parser.set_validate( false );
        parser.parse_memory( documentAsString );
        if (!parser) throw utl::dom::noDocumentParsedXcpt();

        this->attachDocument( parser.get_document() );
    }

    void
    moleculizer::attachDocument(xmlpp::Document* pDoc)
    {
        // Note that this new pattern here will allow old style moleculizer
        // files, ie those that contain events, to be parsed and run.

        if (modelLoaded) throw utl::modelAlreadyLoaded();
        else modelLoaded = true;

        // Get the basic framework of the document.
        xmlpp::Element* pRootElement
            = pDoc->get_root_node();

        // Get the unique model element.
        xmlpp::Element* pModelElement
            = utl::dom::mustGetUniqueChild(pRootElement,
                                           eltName::model);

        // Extract model info.
        constructorCore(pRootElement,
                        pModelElement);

        // Have each unit do its prepareToRun thing.
        std::for_each(pUserUnits->begin(),
                      pUserUnits->end(),
                      prepareUnitToRun(pRootElement,
                                       pModelElement));

    }


    void
    moleculizer::setGenerateDepth(unsigned int generateDepth)
    {
        mzrUnit& rMzrUnit = *(pUserUnits->pMzrUnit);
        rMzrUnit.setGenerateDepth(generateDepth);
    }

     void
     moleculizer::setTolerance(double tolerance)
     {
       // Check that tolerance is a non-negative double.
       mzrReaction::setTolerance( tolerance );
     }

    class unitInsertStateElements :
        public std::unary_function<unit*, void>
    {
        xmlpp::Element* pModelElt;
    public:
        unitInsertStateElements(xmlpp::Element* pModelElement) :
            pModelElt(pModelElement)
        {}

        void
        operator()(unit* pUnit) const throw(std::exception)
        {
            pUnit->insertStateElts(pModelElt);
        }
    };

    // State dump, invoked in a dump-state event.
    xmlpp::Document*
    moleculizer::makeDomOutput(void) throw(std::exception)
    {
        xmlpp::Document* pDoc = new xmlpp::Document();

        // Create the moleculizer-state node.
        xmlpp::Element* pRootElt
            = pDoc->create_root_node(eltName::moleculizerState);

        // Add some mandatory elements.
        xmlpp::Element* pModelElt
            = pRootElt->add_child(eltName::model);
        pModelElt->add_child(eltName::unitsStates);
        pModelElt->add_child(eltName::explicitSpeciesTags);
        pModelElt->add_child(eltName::taggedSpecies);
        pModelElt->add_child(eltName::tagReactions);

        // The only reason for inserting these elements here seems
        // to have been (thought to be?) that various modules would insert
        // elements under them.
        pModelElt->add_child(eltName::time);
        pRootElt->add_child(eltName::streams);

        // Run through the units, letting each make its complete state
        // contribution, for now.
        std::for_each(pUserUnits->begin(),
                      pUserUnits->end(),
                      unitInsertStateElements(pRootElt));

        return pDoc;
    }

    void
    moleculizer::verifyInput(xmlpp::Element const * const pRootElement,
                             xmlpp::Element const * const pModelElement) const
        throw(std::exception)
    {
        // Verify that every element in the model section is handled by some
        // unit or another.
        //
        // I discover very late in the game that this code incorrectly tries
        // to convert comments to elements.
        xmlpp::Node::NodeList modelContentNodes
            = pModelElement->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledModelContent
            = std::find_if(modelContentNodes.begin(),
                           modelContentNodes.end(),
                           modelNodeNotInCap(inputCap));
        if(modelContentNodes.end() != iUnhandledModelContent)
            throw unhandledModelContentXcpt(*iUnhandledModelContent);

        // Get the reaction-gens node, which contains kinds of reaction generators
        // introduced by units.
        xmlpp::Element* pReactionGensElt
            = utl::dom::mustGetUniqueChild(pModelElement,
                                           eltName::reactionGens);

        // Verify that every kind of reaction generator is handled by some
        // unit or another.
        xmlpp::Node::NodeList reactionGensContentNodes
            = pReactionGensElt->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledReactionGenContent
            = std::find_if(reactionGensContentNodes.begin(),
                           reactionGensContentNodes.end(),
                           reactionGenNotInCap(inputCap));
        if(reactionGensContentNodes.end() !=
           iUnhandledReactionGenContent)
            throw unhandledReactionGenXcpt(*iUnhandledReactionGenContent);

        // Get the explicit species model node, which can contain kinds of explicit
        // species introduced by units.
        xmlpp::Element* pExplicitSpeciesElt
            = utl::dom::mustGetUniqueChild(pModelElement,
                                           eltName::explicitSpecies);

        // Verify that every kind of explicit species is handled by some
        // unit or another.
        xmlpp::Node::NodeList explicitSpeciesContentNodes
            = pExplicitSpeciesElt->get_children();
        xmlpp::Node::NodeList::iterator iUnhandledExplicitSpeciesContent
            = std::find_if(explicitSpeciesContentNodes.begin(),
                           explicitSpeciesContentNodes.end(),
                           explicitSpeciesNodeNotInCap(inputCap));

        if(explicitSpeciesContentNodes.end() != iUnhandledExplicitSpeciesContent)
            throw unhandledExplicitSpeciesContentXcpt(*iUnhandledExplicitSpeciesContent);

    }

}

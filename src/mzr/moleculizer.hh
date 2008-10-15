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


#ifndef MOLECULIZER_H
#define MOLECULIZER_H

/*! \defgroup unitsGroup Moleculizer units.
\brief Units and other libraries. */

/*! \defgroup mzrGroup The mzr library.
\ingroup unitsGroup
\brief Library against which main program and all units are linked.

This library is not a moleculizer unit in the same sense as the rest;
it cannot be loaded with the load command (loadCmd).  Instead, it contains
base code for all moleculizer operations, and it is loaded automatically
by the dynamic linker, ld.so.

By linking the main executable and all units against this library,
their sizes are reduced considerably.  If not linked against this
library, the unit .so's are much larger, presumably due to g++'s
current (default) template expansion strategy.  I'm not absolutely
sure that this bloat would correspondingly increase memory
consumption when the libraries are dlopen'ed, but the usual first
step in dlopen is to memory map the .so file, I think. */

/*! \file moleculizer.hh
\ingroup mzrGroup
\brief Defines the main application class. */

/*! \mainpage Moleculizer source code

Originally, Moleculizer units were called "modules," but the
name was changed to avoid conflict with Doxygen's notion of a
module,  a topic-oriented chunk of hierarchical documentation.
You can see a hierarchical view of all the modules by using the
"Modules" item in the page-top menu.  Most of the modules are
in fact Moleculizer units.

Quick links:
- <A HREF="../../index.html">Up to main index.</A>
- \link unitsGroup Moleculizer units. \endlink
*/

#include "utl/defs.hh"
#include "fnd/reactionNetworkCatalog.hh"
#include "mzr/mzrException.hh"
#include "mzr/unit.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrReaction.hh"
#include <boost/foreach.hpp>

namespace mzr
{
    class unitsMgr;

    /*! \ingroup mzrGroup
    \brief The main application object. */


    // The main bulk of this class can be found in ReactionNetworkDescription.
    class moleculizer :
                public fnd::ReactionNetworkDescription<mzrSpecies, mzrReaction>
    {
    public:
        void enableSpatialExtrapolation();
        void enableNonspatialExtrapolation();

    public:
        moleculizer (void);
        ~moleculizer (void);

        bool getModelHasBeenLoaded() const;
        void setModelHasBeenLoaded (bool value);

        void setGenerateDepth (unsigned int generateDepth);
        void setRateExtrapolation ( bool rateExtrapolation )
        {
            return;
        }

        void attachFileName (const std::string& aFileName);
        void attachString (const std::string& documentAsString);
        void attachDocument (xmlpp::Document* pDoc);

        mzrSpecies*
        getSpeciesWithName (const std::string& speciesName) throw ( mzr::IllegalNameXcpt );

    public:
        void
        incrementSpecies (std::string& speciesName);

        std::string
        getRandomSpeciesName() const;

    public:
        xmlpp::Document*
        makeDomOutput (void) throw (std::exception);

    protected:
        void
        constructorPrelude (void);

        void
        verifyInput (const xmlpp::Element* const pRootElt,
                     const xmlpp::Element* const pModelElt,
                     const xmlpp::Element* const pStreamElt) const
            throw (std::exception);

        void
        constructorCore (xmlpp::Element* pRootElt,
                         xmlpp::Element* pModelElt,
                         xmlpp::Element* pStreamElt)
            throw (std::exception);



        // NEW STUFF --------------

        std::map<std::string, std::string> userNameToSpeciesIDChart;

    public: 

        void recordPlexParameter( const std::string& plexName, 
                                  const std::string& parameterName, 
                                  const double& parameterValue);

        void recordUserNameToGeneratedNamePair( const std::string& userName,
                                                const std::string& genName) 
        {
            userNameToSpeciesIDChart.insert( std::make_pair( userName, genName ) );
        }
        
    protected:
        std::string getGeneratedNameFromUserName(const std::string& userName)
        {
            return userNameToSpeciesIDChart[ userName ];
        }

        std::map<SpeciesID, Real> radiusChart; // These are in meters. 
        std::map<SpeciesID, Real> k_DChart;

        std::map<const mzrReaction*, Real> reactionRateChart;
        std::map<const mzrReaction*, Real> binaryActivationEnergyChart;

        std::map<const fnd::coreRxnGen*, double> binaryActivationEnergiesParameterLookup;
        std::map<const fnd::coreRxnGen*, double> unaryReactionRatesParameterLookup;
        std::map<const fnd::coreRxnGen*, double> binaryReactionRatesParameterLookup;

        void installKForNewUnaryReaction( const mzrReaction* mzrReaction);
        void installKForNewBinaryReaction( const mzrReaction* mzrReaction);
        void installKaForNewBinaryReaction( const mzrReaction* mzrReaction);
        void installRadiusForNewSpecies( const mzrReaction* mzrRxn );
        void installKdForNewSpecies( const mzrReaction* mzrRxn );
        double calculateNewRadiusForSpecies( const mzrSpecies* ptrSpecies );

        bool isDecompositionRxn( const mzrReaction* ptrRxn) const
        {
            if (ptrRxn->getNumberOfReactants() == 1 && ptrRxn->getNumberOfProducts() == 2)
                return true;
            else return false;
        }

        bool isDimerizationRxn( const mzrReaction* ptrRxn) const
        {
            if (ptrRxn->getNumberOfReactants() == 2 && ptrRxn->getNumberOfProducts() == 1)
                return true;
            else return false;
        }

        bool isOneToOneRxn( const mzrReaction* ptrRxn) const
        {
            if (ptrRxn->getNumberOfReactants() == 1 && ptrRxn->getNumberOfProducts() == 1)
                return true;
            else 
                return false;
        }

        ////////////////////////////////////////////////
        
    public:
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

    private:
        bool modelLoaded;
    };
}

#endif

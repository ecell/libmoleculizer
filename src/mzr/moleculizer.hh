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

#include "fnd/reactionNetworkCatalog.hh"
#include "utl/autoCatalog.hh"
#include "mzr/unit.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrReaction.hh"
#include "mzr/illegalSpeciesNameXcpt.hh"

namespace mzr
{
    class unitsMgr;
  
    /*! \ingroup mzrGroup
      \brief The main application object. */


    //  The main bulk of this class can be found in ReactionNetworkDescription.
    class moleculizer :
        public fnd::ReactionNetworkDescription<mzrSpecies, mzrReaction>
    {
    public:

        void greeting();

        void RunInteractiveDebugMode();
        void RunProfileMode(unsigned int numIters = 100, bool verbose = false);

        void attachFileName(const std::string& aFileName);
        void attachString(const std::string& documentAsString);
        void attachDocument(xmlpp::Document* pDoc);
      
    public:
       
    public:

        void setGenerateDepth(unsigned int generateDepth);
        void setRateExtrapolation( bool rateExtrapolation ){ return; }
        void setToleranceOption( bool tolerenceOption ) { return; }
        void setTolerance(double tolerance);

    public: 
        // DEBUG interface.  Most of these should be const,
        // but since the InteractiveModeManager stores pointers to member 
        // functions, it gets mad if some are const and others non-const.

        // This is fixable though, with an overloaded addFunction method
        // in InteractiveModeManager and an extra vector of function ptr 
        // storage.

        void 
        incrementSpecies(std::string& speciesName);

        void DEBUG_showTotalNumberSpecies();
        void DEBUG_showTotalNumberReactions();
        void DEBUG_showNumberDeltaSpecies();
        void DEBUG_showNumberDeltaReactions();

        void DEBUG_showAllSpecies();
        void DEBUG_showAllReactions();

        void DEBUG_showDeltaSpecies();
        void DEBUG_showDeltaReactions();

        void DEBUG_incrementSpecies();

        std::string DEBUG_getRandomLiveSpecies() const;
        void DEBUG_incrementRandomSpecies();
        void DEBUG_incrementNRandomLiveSpecies(unsigned int number = 1 );

        void DEBUG_showRandomLiveSpecies();
        void DEBUG_showNumberLiveSpecies();
        void DEBUG_showUnaryReaction();
        void DEBUG_showReactionsBetweenSpecies();
        void DEBUG_changeNamingStrategy();
        void DEBUG_getSpeciesFromName();

        void DEBUG_clearAll();
        void DEBUG_findReaction();
        
        public:
        moleculizer(void);
        ~moleculizer(void);

      
        xmlpp::Document*
        makeDomOutput(void) throw(std::exception);


    protected:
        void
        constructorPrelude(void);

        void 
        verifyInput(xmlpp::Element const * const pRootElt,
                    xmlpp::Element const * const pModelElt,
                    xmlpp::Element const * const pStreamsElt) const 
            throw(std::exception);
    
        void
        constructorCore(xmlpp::Element* pRootElt,
                        xmlpp::Element* pModelElt,
                        xmlpp::Element* pStreamsElt)
            throw(std::exception);


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

    private:
        bool modelLoaded;
    };
}

#endif

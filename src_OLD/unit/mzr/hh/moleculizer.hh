/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001  Walter Lawrence (Larry) Lok.
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

#include "mzr/eventQueue.hh"
#include "mzr/unit.hh"
#include "utl/domJob.hh"

namespace mzr
{
  class unitsMgr;
  
  /*! \ingroup mzrGroup
    \brief The main application object. */
  class moleculizer :
    public utl::dom::domBatchApp
  {
  protected:
    void
    constructorPrelude(void);
    
    void
    constructorCore(xmlpp::Element* pRootElt,
		    xmlpp::Element* pModelElt,
		    xmlpp::Element* pStreamsElt,
		    xmlpp::Element* pEventsElt)
      throw(std::exception);

  public:

    // Units loaded by the user, waiting for destruction.
    // 
    // This class is the manager for units, and the place that new units
    // can be installed.  It's public because units need to get to
    // each other.
    unitsMgr* pUserUnits;

    eventQueue eventQ;

    // Codes the input capabilities of moleculizer, including its
    // parsing routine.
    inputCapabilities inputCap;

    // This constructor is used in the construction of the two
    // descendant apps, continuator and parametrizer.
    moleculizer(void)
      throw(std::exception);

    // This constructor is used for runs as moleculizer.
    moleculizer(int argc,
		char** argv,
		xmlpp::Document* pDoc) throw(std::exception);

    ~moleculizer(void);

    void
    processCommandLineArgs(int argc,
			   char* argv[]);
    
    // Were it not for generating output with xmlpp code, which throws base
    // std::exceptions, I would probably just be throwing my own exceptions
    // here.
    int
    run(void) throw(std::exception);

    // Generate an output document for rk4tau or rk4ode.
    xmlpp::Document*
    makeDomOutput(void) throw(std::exception);
  };
}

#endif

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

#ifndef CPT_CPTUNIT_H
#define CPT_CPTUNIT_H

#include "utl/dom.hh"
#include "utl/gsl.hh"
#include "utl/autoCatalog.hh"
#include "utl/autoVector.hh"

#include "fnd/query.hh"

#include "cpt/unit.hh"
#include "cpt/tabDumpEvent.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/singleGlobalSpeciesDumpable.hh"
#include "cpt/dumpableNotSpeciesStreamXcpt.hh"

namespace cpt
{
  class cptUnit :
    public unit
  {
    // The overall structure of space.  All of the compartments' volumes
    // are state variables, which lie down inside this thing.
    compartmentGraph theGraph;

    // List/vector/??? of all global species, for traversal during diffusion.
    std::vector<globalSpecies*> allGlobalSpecies;

    // These are all the species that can be dumped to output.
    utl::catalog<globalSpecies> speciesByName;

    // All the species that are not deleted by the recognizer and
    // therefore need management here.
    utl::autoVector<globalSpecies> managedSpecies;

    // Memory management of reaction families, which hold all the automatically
    // generated reactions.
    //
    // Are reaction families really of any use?  A reaction family consists of
    // all the reactions made by a particular reaction generator, and formerly
    // this classification of reactions was preserved in state dump output.
    // That is, you could tell which reactions in a state dump were generated
    // by which reaction generators.
    utl::autoVector<std::vector<globalReaction*> > reactionFamilies;

    // I see below that these require special initial scheduling???
    // These are also in userEvents for memory management.
    std::vector<tabDumpEvent*> tabDumpEvents;

    // All dumpables.
    utl::autoCatalog<fnd::basicDumpable>dumpables;

    // Each entry here should correspond to an entry in the dumpables
    // catalog with the same name.
    std::vector<cptSpeciesStream*> speciesStreams;

    // Entries here are mainly mixins in tabDumpEvents, and represent .dmp
    // files.  With the speciesStreams this arranges that one can connect
    // .dmp output of populations with species in a state dump.  The keys
    // here are basically output file name strings.
    std::vector<cptDumpStream*> dumpStreams;

    // Memory management for queries of all kinds.  The basicQuery
    // base class serves no purpose other than memory management here.
    utl::autoVector<fnd::baseQuery> queries;

    // The one and only random number generator.  Reseeding this reseeds
    // everything.
    utl::gsl::autoGslRng rng;

    // These are connected with the "-g" command line options, which
    // turns network generation off.  This allows one to load a simulation
    // state, then simulate as an ordinary stochastic simulator.
    //
    // To do this, one needs have automatic generation on, when the
    // simulation state is loaded, and then turned off if indicated by
    // the command line option.  See msrUnit::prepareToRun.
    bool generateOption;
    bool generateOk;

    // Connected with the "-d" option, this is the depth to which species
    // notify in the reaction network.  This allows one to do "generate, then
    // simulate" mode with Moleculizer, for example.  It is perhaps more
    // intereseting to just make a Moleculizer-generated network a little
    // "thicker" by setting this to a strictly positive value.
    int generateDepth;

  public:

    cptUnit(cptApp& rCptApp);

    compartmentGraph&
    getCompartmentGraph(void)
    {
      return theGraph;
    }

    const compartmentGraph&
    getCompartmentGraph(void) const
    {
      return theGraph;
    }

    // The time this moleculizer was created; taken to be the start
    // time of the program.  Used in secondsDumpable.  This being static
    // seems ridiculous; it's easier to get to in the dumpable code, though.
    time_t startSeconds;

    // Connected with the startSeconds, this member variable can limit
    // moleculizer to run only a certain number of clock seconds.  If
    // this number is 0, no check is performed.  If positive, clock seconds
    // are checked against it on every event(!) and event processing
    // is terminated when the time is up.
    time_t secondsLimit;

    void
    setTimeout(time_t seconds)
    {
      secondsLimit = seconds;
    }

    // Adds the the species to the species that are tranversed during
    // diffusion.  This should effectively be all global species.
    //
    // It might make sense to put registration here in the globalSpecies
    // constructor, but if that is done, then the other "addXXXSpecies"
    // routines below will need to be changed.
    void
    addSpecies(globalSpecies* pSpecies)
    {
      allGlobalSpecies.push_back(pSpecies);
    }

    // Access to all global species for diffusion.
    const std::vector<globalSpecies*>&
    getGlobalSpecies(void) const
    {
      return allGlobalSpecies;
    }

    // Adds a species that has a name, so that it is has a corresponding
    // dumpable.  Returns false if the name duplicates the name of another
    // species or another dumpable.
    bool
    addNamedSpecies(const std::string& rSpeciesName,
		    globalSpecies* pSpecies)
    {
      bool nameOk 
	= speciesByName.addEntry(rSpeciesName,
				 pSpecies)
	&& addSpeciesDumpable(new singleGlobalSpeciesDumpable(rSpeciesName,
							      pSpecies));
      if(nameOk) addSpecies(pSpecies);
      
      return nameOk;
    }

    // Adds the named species, throwing an exception if the name
    // duplicates the name of another species or another dumpable.
    void
    mustAddNamedSpecies(const std::string& rSpeciesName,
			globalSpecies* pSpecies,
			xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // Gives a name to a species that has already been added.
    // This is a stupid hack to get around the fact that plexSpecies
    // are "added" when they are created by their plexFamily, but there
    // are named ("explicit") plexSpecies.
    void
    mustNameSpecies(const std::string& rSpeciesName,
		    globalSpecies* pSpecies,
		    xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    globalSpecies*
    findSpecies(const std::string& rSpeciesName) const
    {
      return speciesByName.findEntry(rSpeciesName);
    }

    globalSpecies*
    mustFindSpecies(const std::string& rSpeciesName,
		    xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);


    // For committing a species to memory management by cptUnit.
    //
    // No longer in use?
    void
    manageSpecies(globalSpecies* pSpecies)
    {
      managedSpecies.push_back(pSpecies);
    }

    // Memory management and traversal of reactions that are not automatically
    // generated.
    utl::autoVector<globalReaction> userReactions;

    void
    addUserReaction(globalReaction* pReaction)
    {
      userReactions.push_back(pReaction);
    }

    // Memory management and traversal of reactions that are automatically
    // generated.
    //
    // It occurs to me that the only reason for these to exist is that the
    // reactions in a family are generated by a particular reaction generator.
    // It might therefore make sense to include in reactionFamily a reference
    // back to its generator.  One could even link generators back to user
    // input for diagnostics, etc.  Formerly, reaction family info was
    // included in state dump.
    void
    addReactionFamily(utl::autoVector<globalReaction>* pReactionFamily)
    {
      reactionFamilies.push_back(pReactionFamily);
    }
  
    // Memory management of events created by the user, such as the stop event.
    // Reactions aren't even events in this sense in cpt.
    utl::autoVector<cptEvent> userEvents;

    void
    addUserEvent(cptEvent* pEvent)
    {
      userEvents.push_back(pEvent);
    }

    // Adds a tabDumpEvent both to the special list for state dumping
    // and for memory management.
    void
    addTabDumpEvent(tabDumpEvent* pTabDumpEvent)
    {
      addUserEvent(pTabDumpEvent);
      
      tabDumpEvents.push_back(pTabDumpEvent);
    }

    bool
    addDumpable(fnd::basicDumpable* pDumpable)
    {
      return dumpables.addEntry(pDumpable->getName(),
				pDumpable);
    }

    // Throws an exception if there is already a dumpable
    // whose name duplicates that of the given dumpable.
    void
    mustAddDumpable(fnd::basicDumpable* pDumpable,
		    xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // The dumpables that print species populations (speciesStreams) are
    // tracked, so that when state is dumped, those species can be enumerated.
    template<class dumpableType>
    bool
    addSpeciesDumpable(dumpableType* pDumpable)
    {
      cptSpeciesStream* pSpeciesStream
	= dynamic_cast<cptSpeciesStream*>(pDumpable);

      if(! pSpeciesStream)
	throw dumpableNotSpeciesStreamXcpt(pDumpable->getName());

      bool nameOk = addDumpable(pDumpable);
      if(nameOk) addSpeciesStream(pSpeciesStream);

      return nameOk;
    }

    fnd::basicDumpable*
    findDumpable(const std::string& rDumpableName) const
    {
      return dumpables.findEntry(rDumpableName);
    }

    fnd::basicDumpable*
    mustFindDumpable(const std::string& rDumpableName,
		     xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    void
    addSpeciesStream(cptSpeciesStream* pStream)
    {
      speciesStreams.push_back(pStream);
    }
    
    // A dumpStream represents a .dmp output file.
    // Here the "streamName" is the file name, +, or -.
    void
    addDumpStream(cptDumpStream* pStream)
    {
      dumpStreams.push_back(pStream);
    }

    void
    addQuery(fnd::baseQuery* pQuery)
    {
      queries.push_back(pQuery);
    }

    // Determine if species/reaction generation is on or off.
    bool
    getGenerateOk(void) const
    {
      return generateOk;
    }

    // Turn reaction generation off or on.  The "noReactEvent" turns (turned?)
    // reaction generation off, but doesn't currently seem to have any input
    // syntax.
    void
    setGenerateOk(bool newGenerateOk)
    {
      generateOk = newGenerateOk;
    }

    // Returns the depth to probe in reaction network generation.
    int
    getGenerateDepth(void)
    {
      return generateDepth;
    }

    // Note that this is pointlessly exposing a non-constant reference.
    // rng could just as well be a public member variable.
    utl::gsl::autoGslRng&
    getRng(void)
    {
      return rng;
    }

    // As in Moleculizer, I'm making this library the home of a number
    // of globals so that they can be shared by different applications
    // (as with moleculizer, continuator, parametrizer).
    //
    // Now it seems to me that there should just be one executable with
    // different functions selected by command line argument.
    //
    // At any rate, for now the command line arguments are such globals.
    void
    parseCommandLine(int argc,
		     char* argv[]);

    virtual void
    parseDomInput(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt,
		  xmlpp::Element* pEventsElt) throw(std::exception);

    // Just emits header lines in all the dumpables, schedules tabDumpEvents
    // for the first time, after which they schedule themselves.
    void
    prepareToRun(xmlpp::Element* pRootElt,
		 xmlpp::Element* pModelElt,
		 xmlpp::Element* pStreamsElt,
		 xmlpp::Element* pEventsElt) throw(std::exception);

    void
    prepareToDump(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt,
		  xmlpp::Element* pEventsElt,
		  xmlpp::Element* pTaggedSpeciesElement)
      throw(std::exception)
    {
      // The default here would be "prepareToRun" which opens dumpfiles and
      // schedules tabDumpEvents for the first time.  This is not called
      // for here.  In fact, the working directory where we're running
      // parametrizer isn't usually writeable.
    }

    // In addition to the above, sets the current (i.e. initial)
    // simulation time to the time at which state was dumped.
    void
    prepareToContinue(xmlpp::Element* pRootElt,
		      xmlpp::Element* pModelElt,
		      xmlpp::Element* pStreamsElt,
		      xmlpp::Element* pEventsElt,
		      std::map<std::string, std::string>& rTagToName,
		      xmlpp::Element* pTaggedSpeciesElement)
      throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception);
  };
}

#endif // CPT_CPTUNIT_H

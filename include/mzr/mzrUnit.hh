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

#ifndef MZRUNIT_H
#define MZRUNIT_H

#include "mzr/util.hh"
#include "domUtils/domUtils.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/reactionFamily.hh"
#include "mzr/unit.hh"
#include "mzr/molarFactor.hh"
#include "mzr/event.hh"
#include "mzr/species.hh"
#include "mzr/mzrUnitParse.hh"

namespace mzr
{
  class unitsMgr;

  // This unit is vaguely special in that, as a shared object, it contains
  // the code for the application class, moleculizer.  As a parser,
  // it "cleans up" general stuff that is basic to moleculizer's operation,
  // rather than added in other units.
  class mzrUnit : public unit
  {
    // The one (?) built in state variable, essentially the volume.
    molarFactor theMolarFactor;

    // These are all the species that can be dumped to output.
    catalog<species> speciesByName;

    // All the species that are not deleted by the recognizer and
    // therefore need management here.
    //
    // 29May03 I will probably want to farm these species out to
    // their units.  For the time being, these are probably just stochastirator
    // species.  Any chance of just keeping keeping, say stochastirator
    // species in a vector, to automate their deletion more gracefully?
    //
    // Each different kind of species would need its own list, and each
    // different kind of species would have to have a default constructor
    // for this kind of thing to work.
    autoVector<species> userSpecies;

    // Memory management of reaction families, which hold all the automatically
    // generated reactions.
    autoVector<reactionFamily> reactionFamilies;

    // Events that were created by the user and therefore need deletion.
    autoVector<event> userEvents;

    // Tab-dump events, which play a special role in state dump.
    // These are also in userEvents for memory management.
    std::vector<tabDumpEvent*> tabDumpEvents;

    // All dumpables.
    autoCatalog<dumpable> userDumpables;

    // Species dumpables.  These also play a special role in state dump;
    // only speciesDumpables are dumped as "speciesStreams".
    catalog<speciesDumpable> speciesDumpables;
    
  public:

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

    // These can be used to turn off the generation of new reactions, and
    // therefore new species.  Each reaction generator looks at generateOk
    // when notified of a new species before generating reactions
    // for that new species.
    //
    // For now, noGenerateOption is set when the program starts up, and copied
    // to generateOk when simulation starts.  This allows the reactions
    // associated with explicit species to be generated when the program
    // starts, but no further reactions.  See mzrUnit::prepareToRun.
    bool generateOption;
    bool generateOk;

    molarFactor&
    getMolarFactor(void)
    {
      return theMolarFactor;
    }

    // This doesn't register the species for deletion, so it can be used
    // for the "summary species" of a plexFamily.
    bool
    addSpecies(const std::string& rSpeciesName,
	       species* pSpecies)
    {
      return speciesByName.addEntry(rSpeciesName,
				    pSpecies)
	&& userDumpables.addEntry(rSpeciesName,
				  new singleSpeciesDumpable(rSpeciesName,
							    *pSpecies));
    }

    species*
    findSpecies(const std::string& rSpeciesName)
    {
      return speciesByName.findEntry(rSpeciesName);
    }

    species*
    mustFindSpecies(xmlpp::Node* pRequestingNode,
		    const std::string& rSpeciesName)
      throw(unknownSpeciesXcpt);

    // This is for adding a species that will be deleted by userData,
    // This is NOT for species that will take care of their own deletion,
    // such as the "summary" species of a plexFamily.
    bool
    addUserSpecies(const std::string& rSpeciesName,
		   species* pSpecies)
    {
      bool result = addSpecies(rSpeciesName,
			       pSpecies);

      if(result) userSpecies.push_back(pSpecies);

      return result;
    }

    // This is for data recursion to do full Levchenko-style
    // extrapolation of the reaction network.  It is not connected with
    // memory management.
    std::vector<reaction*>allReactions;

    void
    addReaction(reaction* pReaction)
    {
      allReactions.push_back(pReaction);
    }

    // Memory management of reactions that are not automatically
    // generated.
    autoVector<reaction> userReactions;

    void
    addUserReaction(reaction* pReaction)
    {
      userReactions.push_back(pReaction);
      addReaction(pReaction);
    }

    void
    addReactionFamily(reactionFamily* pReactionFamily)
    {
      reactionFamilies.push_back(pReactionFamily);
    }
  
    void
    addEvent(event* pEvent)
    {
      userEvents.push_back(pEvent);
    }

    // These events play a special role in state dump.
    void
    addTabDumpEvent(tabDumpEvent* pEvent)
    {
      tabDumpEvents.push_back(pEvent);
      addEvent(pEvent);
    }

    bool
    addDumpable(dumpable* pDumpable)
    {
      return userDumpables.addEntry(pDumpable->getName(),
				    pDumpable);
    }

    // These play a special role in state dump.
    bool
    addSpeciesDumpable(speciesDumpable* pSpeciesDumpable)
    {
      return addDumpable(pSpeciesDumpable)
	&& speciesDumpables.addEntry(pSpeciesDumpable->getName(),
				     pSpeciesDumpable);
    }

    dumpable*
    findDumpable(const std::string& rDumpableName)
    {
      return userDumpables.findEntry(rDumpableName);
    }

    speciesDumpable*
    findSpeciesDumpable(const std::string& rDumpableName)
    {
      return speciesDumpables.findEntry(rDumpableName);
    }

    mzrUnit(moleculizer& rMoleculizer) :
      unit("mzr",
	   rMoleculizer),
      startSeconds(time(0)),
      secondsLimit(0),
      generateOption(true),
      generateOk(true)
    {
    
      // Model elements whose contents are parsed by some unit
      // or another, as determined by moleculizer::parseDomInput.
      inputCap.addModelContentName(eltName::reactionGens);
      inputCap.addModelContentName(eltName::explicitSpecies);

      // Model elements that this unit actually processes.
      inputCap.addModelContentName(eltName::explicitReactions);
      inputCap.addModelContentName(eltName::volume);

      // This unit is not responsible for any reaction generators
      // or species streams.

      // Events for which this unit is responsible.
      inputCap.addEventsContentName(eltName::createEvent);
      inputCap.addEventsContentName(eltName::dumpStateEvent);
      inputCap.addEventsContentName(eltName::stopEvent);

      // Add the standard, built-in dumpables.
      addDumpable(new clockDumpable());
      addDumpable(new secondsDumpable(*this));
      addDumpable(new reactEventCountDumpable());
      addDumpable(new reactCountDumpable());
      addDumpable(new activationCountDumpable());
      addDumpable(new speciesCountDumpable());
      addDumpable(new simTimeDumpable(rMolzer));
      addDumpable(new volumeDumpable(*this));
    }

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
    mzrUnit::prepareToDump(xmlpp::Element* pRootElt,
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

#endif // MZRUNIT_H

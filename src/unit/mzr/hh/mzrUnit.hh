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

#include "utl/dom.hh"
#include "utl/gsl.hh"
#include "utl/autoCatalog.hh"
#include "utl/autoVector.hh"

#include "mzr/unit.hh"
#include "mzr/molarFactor.hh"
#include "mzr/mzrEvent.hh"
#include "mzr/tabDumpEvent.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrSpeciesDumpable.hh"
#include "mzr/dumpableNotSpeciesStreamXcpt.hh"

namespace mzr
{
  // This unit is vaguely special in that, as a shared object, it contains
  // the code for the application class, moleculizer.  As a parser,
  // it "cleans up" general stuff that is basic to moleculizer's operation,
  // rather than added in other units.
  class mzrUnit :
    public unit
  {
    // The one (?) built in state variable, essentially the volume.
    molarFactorGlobal theMolarFactor;

    // These are all the species that can be dumped to output.
    utl::catalog<mzrSpecies> speciesByName;

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
    utl::autoVector<mzrSpecies> userSpecies;

    // Memory management of reaction families, which hold all the automatically
    // generated reactions.
    utl::autoVector<utl::autoVector<mzrReaction> > reactionFamilies;

    // I see below that these require special initial scheduling???
    // These are also in userEvents for memory management.
    std::vector<tabDumpEvent*> tabDumpEvents;

    // All dumpables.
    utl::autoCatalog<fnd::dumpable<fnd::basicDumpable::dumpArg> >dumpables;

    // Each entry here should correspond to an entry in the dumpables
    // catalog with the same name.
    std::vector<mzrSpeciesStream*> speciesStreams;

    // Entries here are mainly mixins in tabDumpEvents, and represent .dmp
    // files.  With the speciesStreams this arranges that one can connect
    // .dmp output of populations with species in a state dump.  The keys
    // here are basically output file name strings.
    std::vector<mzrDumpStream*> dumpStreams;

    // Memory management for queries of all kinds.  The basicQuery
    // base class serves no purpose other than memory management here.
    utl::autoVector<fnd::baseQuery> queries;

    // This new command-line option supplants generateOption and generateOk,
    // in that making the depth negative should turn off reaction generation
    // after intial setup.
    int generateDepth;

    bool generateOption;
    bool generateOk;

  public:

    // Iterators into this are used in constructor of reaction to ensure that
    // each reaction is sensitized to each global state variable.
    std::vector<fnd::sensitivityList<mzrReaction>*> globalVars;

    mzrUnit(moleculizer& rMoleculizer);

    // The one and only random number generator.
    utl::gsl::autoGslRng rng;

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

    // Accessors for generation depth command-line argument.
    void
    setGenerateDepth(int depth)
    {
      generateDepth = depth;
    }

    int
    getGenerateDepth(void)
    {
      return generateDepth;
    }

    bool
    getGenerateOk(void)
    {
      return generateOk;
    }

    void
    setGenerateOk(bool isGenerateOk)
    {
      generateOk = isGenerateOk;
    }

    bool
    getGenerateOption(void)
    {
      return generateOption;
    }

    void
    setGenerateOption(bool isGenerationOn)
    {
      generateOption = isGenerationOn;
    }

    molarFactorGlobal&
    getMolarFactor(void)
    {
      return theMolarFactor;
    }

    // This doesn't register the species for deletion, so it can be used
    // for explicit plexSpecies, which are deleted by their plexFamilies.
    bool
    addSpecies(const std::string& rSpeciesName,
	       mzrSpecies* pSpecies)
    {
      return speciesByName.addEntry(rSpeciesName, (mzrSpecies*) pSpecies)
              && dumpables.addEntry(rSpeciesName,
	      new singleSpeciesDumpable<mzrSpecies>(rSpeciesName,pSpecies));
    }

    void
    mustAddSpecies(const std::string& rSpeciesName,
		   mzrSpecies* pSpecies,
		   xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // For adding a species that will be memory-managed by mzrUnit.
    // (An example would be an explicit stochSpecies, but that's from
    // another module.)
    bool
    addUserSpecies(const std::string& rSpeciesName,
		   mzrSpecies* pSpecies)
    {
      bool result = addSpecies(rSpeciesName,
			       pSpecies);

      if(result) userSpecies.push_back(pSpecies);

      return result;
    }

    mzrSpecies*
    findSpecies(const std::string& rSpeciesName) const
    {
      return speciesByName.findEntry(rSpeciesName);
    }

    mzrSpecies*
    mustFindSpecies(const std::string& rSpeciesName,
		    xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);


    // Memory management and traversal of reactions that are not automatically
    // generated.
    utl::autoVector<mzrReaction> userReactions;

    void
    addUserReaction(mzrReaction* pReaction)
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
    // input for diagnostics, etc.
    void
    addReactionFamily(utl::autoVector<mzrReaction>* pReactionFamily)
    {
      reactionFamilies.addEntry(pReactionFamily);
    }
  
    // Memory management of events created by the user.
    utl::autoVector<mzrEvent> userEvents;

    void
    addUserEvent(mzrEvent* pEvent)
    {
      userEvents.push_back(pEvent);
    }

    bool
    addDumpable(fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable)
    {
      return dumpables.addEntry(pDumpable->getName(),
				pDumpable);
    }

    // Throws an exception if there is already a dumpable
    // whose name duplicates that of the given dumpable.
    void
    mustAddDumpable(fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable,
		    xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // The dumpables that print species populations (speciesStreams) are
    // tracked, so that when state is dumped, those species can be enumerated.
    template<class dumpableType>
    bool
    addSpeciesDumpable(dumpableType* pDumpable)
    {
      mzrSpeciesStream* pSpeciesStream
	= dynamic_cast<mzrSpeciesStream*>(pDumpable);

      if(! pSpeciesStream)
	throw dumpableNotSpeciesStreamXcpt(pDumpable->getName());

      bool nameOk = addDumpable(pDumpable);
      if(nameOk) addSpeciesStream(pSpeciesStream);

      return nameOk;
    }

    fnd::dumpable<fnd::basicDumpable::dumpArg>*
    findDumpable(const std::string& rDumpableName) const
    {
      return dumpables.findEntry(rDumpableName);
    }

    fnd::dumpable<fnd::basicDumpable::dumpArg>*
    mustFindDumpable(const std::string& rDumpableName,
		     xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    void
    addSpeciesStream(mzrSpeciesStream* pStream)
    {
      speciesStreams.push_back(pStream);
    }
    
//     void
//     mustAddSpeciesStream(xmlpp::Node* pRequestingNode,
// 			 const std::string& rStreamName,
// 			 mzrSpeciesStream* pStream)
//       throw(utl::xcpt);

//     mzrSpeciesStream*
//     findSpeciesStream(const std::string& rStreamName)
//     {
//       return speciesStreams.findEntry(rStreamName);
//     }

//     mzrSpeciesStream*
//     mustFindSpeciesStream(xmlpp::Node* pRequestingNode,
// 			  const std::string& rStreamName)
//       throw(utl::xcpt);

    // A dumpStream represents a .dmp output file.
    // Here the "streamName" is the file name, +, or -.
    void
    addDumpStream(mzrDumpStream* pStream)
    {
      dumpStreams.push_back(pStream);
    }

//     void
//     mustAddDumpStream(const std::string& rStreamName,
// 		      mzrDumpStream* pStream,
// 		      xmlpp::Node* pRequestingNode = 0)
//       throw(utl::xcpt);

//     mzrDumpStream*
//     findDumpStream(const std::string& rStreamName)
//     {
//       return dumpStreams.findEntry(rStreamName);
//     }

//     mzrDumpStream*
//     mustFindDumpStream(const std::string& rStreamName,
// 		       xmlpp::Node* pRequestingNode = 0)
//       throw(utl::xcpt);

//     mzrUnit(moleculizer& rMoleculizer);

    void
    addQuery(fnd::baseQuery* pQuery)
    {
      queries.push_back(pQuery);
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

#endif // MZRUNIT_H

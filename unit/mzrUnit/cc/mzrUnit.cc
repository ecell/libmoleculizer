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

#include "mzr/unitsMgr.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/moleculizer.hh"
#include "mzr/dumpUtils.hh"

namespace mzr
{
  species*
  mzrUnit::
  mustFindSpecies(xmlpp::Node* pRequestingNode,
		  const std::string& rSpeciesName) const
    throw(unknownSpeciesXcpt)
  {
    species* pSpecies = findSpecies(rSpeciesName);

    if(! pSpecies)
      throw unknownSpeciesXcpt(pRequestingNode,
			       rSpeciesName);
    return pSpecies;
  }

  class insertReaction :
    public std::unary_function<const reaction*, void>
  {
    xmlpp::Element* pTagReactionsElt;
  public:
    insertReaction(xmlpp::Element* pTagReactionsElement) :
      pTagReactionsElt(pTagReactionsElement)
    {
    }
    void
    operator()(const reaction* pReaction) throw(std::exception)
    {
      pReaction->insertElt(pTagReactionsElt);
    }
  };

  class insertFamilyReactions :
    public std::unary_function<const reactionFamily*, void>
  {
    // This may change to a different element containing just the
    // automatically generated reactions.
    xmlpp::Element* pTagReactionsElt;
  public:
    insertFamilyReactions(xmlpp::Element* pTagReactionsElement) :
      pTagReactionsElt(pTagReactionsElement)
    {}

    void
    operator()(const reactionFamily* pFamily) const throw(std::exception)
    {
      std::for_each(pFamily->begin(),
		    pFamily->end(),
		    insertReaction(pTagReactionsElt));
    }
  };

  class insertExplicitSpeciesTag :
    public std::unary_function<catalog<species>::value_type, void>
  {
    xmlpp::Element* pExplicitSpeciesTagsElt;
  public:
    insertExplicitSpeciesTag(xmlpp::Element* pExplicitSpeciesTagsElement) :
      pExplicitSpeciesTagsElt(pExplicitSpeciesTagsElement)
    {}

    void
    operator()(const argument_type& rNameSpeciesPair) const
      throw(std::exception)
    {
      xmlpp::Element* pExplicitSpeciesTagElt
	= pExplicitSpeciesTagsElt->add_child(eltName::explicitSpeciesTag);

      pExplicitSpeciesTagElt
	->set_attribute(eltName::explicitSpeciesTag_nameAttr,
			rNameSpeciesPair.first);

      pExplicitSpeciesTagElt
	->set_attribute(eltName::explicitSpeciesTag_tagAttr,
			rNameSpeciesPair.second->getTag());
    }
  };

  // Inserts a tagged-species-stream element for each dumpable
  // in the autoCatalog that is actually a species dumpable.
  //
  // This is an isolated, non-time-critical, but somewhat dubious use of
  // dynamic_cast, since I could maintain a separate list of the species
  // dumpables.
  class insertTaggedSpeciesStream :
    public std::unary_function<autoCatalog<dumpable>::value_type, void>
  {
    xmlpp::Element* pTaggedSpeciesStreamsElt;
  public:
    insertTaggedSpeciesStream(xmlpp::Element* pTaggedSpeciesStreamsElement) :
      pTaggedSpeciesStreamsElt(pTaggedSpeciesStreamsElement)
    {
    }

    void
    operator()(const argument_type& rNameDumpablePair) const
      throw(std::exception)
    {
      speciesDumpable* pSpeciesDumpable
	= dynamic_cast<speciesDumpable*>(rNameDumpablePair.second);

      if(0 != pSpeciesDumpable)
	{
	  xmlpp::Element* pTaggedSpeciesStreamElt
	    = pTaggedSpeciesStreamsElt
	    ->add_child(eltName::taggedSpeciesStream);

	  pTaggedSpeciesStreamElt
	    ->set_attribute(eltName::taggedSpeciesStream_nameAttr,
			    pSpeciesDumpable->getName());

	  pSpeciesDumpable->insertDumpedSpeciesTags(pTaggedSpeciesStreamElt);
	}
    }
  };

  class insertTaggedDumpStream :
    public std::unary_function<const event*, void>
  {
    xmlpp::Element* pTaggedDumpStreamsElt;
  public:
    insertTaggedDumpStream(xmlpp::Element* pTaggedDumpStreamsElement) :
      pTaggedDumpStreamsElt(pTaggedDumpStreamsElement)
    {
    }

    void
    operator()(const event* pEvent) const
      throw(std::exception)
    {
      // This is another convenient, if dubious use of dynamic_cast,
      // since I could keep a separate vector of tabDumpEvents, rather
      // than filtering for them here.
      const tabDumpEvent* pTabDumpEvent
	= dynamic_cast<const tabDumpEvent*>(pEvent);

      if(0 != pTabDumpEvent)
	{
	  pTabDumpEvent->insertTaggedDumpStreamElts(pTaggedDumpStreamsElt);
	}
    }
  };

  void
  mzrUnit::insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
  {
    // Model elements.
    xmlpp::Element* pModelElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     eltName::model);

    xmlpp::Element* pExplicitSpeciesTagsElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     eltName::explicitSpeciesTags);

    // Give tags for named species.
    std::for_each(speciesByName.begin(),
		  speciesByName.end(),
		  insertExplicitSpeciesTag(pExplicitSpeciesTagsElt));

    // Give all the reactions, using tags to refer to species.
    xmlpp::Element* pTagReactionsElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     eltName::tagReactions);

    std::for_each(userReactions.begin(),
		  userReactions.end(),
		  insertReaction(pTagReactionsElt));

    std::for_each(reactionFamilies.begin(),
		  reactionFamilies.end(),
		  insertFamilyReactions(pTagReactionsElt));

    // Insert the volume.  This has to have separate fraction and exponent
    // for SBML scientific notation.
    double volume = getMolarFactor().getVolume();
    addDoubleParamChild(pRootElt,
			eltName::volume,
			eltName::volume_litersAttr,
			volume);

    // Insert the simulation time (when the dump happens) so that
    // the simulation can be continuted from the dump.
    xmlpp::Element* pTimeElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     eltName::time);

    double time = rMolzer.eventQ.getSimTime();
    pTimeElt->set_attribute(eltName::time_secondsAttr,
			    domUtils::stringify<double>(time));

    // Streams element.
    xmlpp::Element* pStreamsElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     eltName::streams);

    // Run through the dumpables, determining which of them dump species,
    // and have those emit the tags of the species that they dump.
    xmlpp::Element* pTaggedSpeciesStreamsElt
      = pStreamsElt->add_child(eltName::taggedSpeciesStreams);
    std::for_each(userDumpables.begin(),
		  userDumpables.end(),
		  insertTaggedSpeciesStream(pTaggedSpeciesStreamsElt));

    // Run through the events, determining which of them are tabDumpEvents,
    // and have those write themselves as dump-streams.  Dump streams in
    // rk4tau have only one name, the file name, a convention that might work
    // well in Moleculizer with the added assumption that dump streams start
    // at time 0.
    xmlpp::Element* pTaggedDumpStreamsElt
      = pStreamsElt->add_child(eltName::taggedDumpStreams);
    std::for_each(userEvents.begin(),
		  userEvents.end(),
		  insertTaggedDumpStream(pTaggedDumpStreamsElt));
  }
}

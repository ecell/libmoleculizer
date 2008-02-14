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

#include "mzr/mzrUnit.hh"
#include "mzr/dumpUtils.hh"
#include "mzr/moleculizer.hh"

namespace mzr
{
//   class insertReaction :
//     public std::unary_function<const mzrReaction*, void>
//   {
//     xmlpp::Element* pTagReactionsElt;
//   public:
//     insertReaction(xmlpp::Element* pTagReactionsElement) :
//       pTagReactionsElt(pTagReactionsElement)
//     {
//     }
//     void
//     operator()(const mzrReaction* pReaction) throw(std::exception)
//     {
//       pReaction->insertElt(pTagReactionsElt);
//     }
//   };

  class insertFamilyReactions :
    public std::unary_function<const std::vector<mzrReaction>*, void>
  {
    // This may change to a different element containing just the
    // automatically generated reactions.
    xmlpp::Element* pTagReactionsElt;
  public:
    insertFamilyReactions(xmlpp::Element* pTagReactionsElement) :
      pTagReactionsElt(pTagReactionsElement)
    {}

    void
    operator()(const std::vector<mzrReaction*>* pFamily) const 
      throw(std::exception)
    {
      std::for_each(pFamily->begin(),
		    pFamily->end(),
		    std::bind2nd(std::mem_fun(&mzrReaction::insertElt),
				 pTagReactionsElt));
    }
  };

  class insertExplicitSpeciesTag :
    public std::unary_function<utl::catalog<mzrSpecies>::value_type, void>
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
//   class insertTaggedSpeciesStream :
//     public std::unary_function<utl::autoCatalog<dumpable>::value_type, void>
//   {
//     xmlpp::Element* pTaggedSpeciesStreamsElt;
//   public:
//     insertTaggedSpeciesStream(xmlpp::Element* pTaggedSpeciesStreamsElement) :
//       pTaggedSpeciesStreamsElt(pTaggedSpeciesStreamsElement)
//     {
//     }

//     void
//     operator()(const argument_type& rNameDumpablePair) const
//       throw(std::exception)
//     {
//       speciesDumpable* pSpeciesDumpable
// 	= dynamic_cast<speciesDumpable*>(rNameDumpablePair.second);

//       if(0 != pSpeciesDumpable)
// 	{
// 	  xmlpp::Element* pTaggedSpeciesStreamElt
// 	    = pTaggedSpeciesStreamsElt
// 	    ->add_child(eltName::taggedSpeciesStream);

// 	  pTaggedSpeciesStreamElt
// 	    ->set_attribute(eltName::taggedSpeciesStream_nameAttr,
// 			    pSpeciesDumpable->getName());

// 	  pSpeciesDumpable->insertDumpedSpeciesTags(pTaggedSpeciesStreamElt);
// 	}
//     }
//   };

//   class insertTaggedDumpStream :
//     public std::unary_function<const event*, void>
//   {
//     xmlpp::Element* pTaggedDumpStreamsElt;
//   public:
//     insertTaggedDumpStream(xmlpp::Element* pTaggedDumpStreamsElement) :
//       pTaggedDumpStreamsElt(pTaggedDumpStreamsElement)
//     {
//     }

//     void
//     operator()(const event* pEvent) const
//       throw(std::exception)
//     {
//       // This is another convenient, if dubious use of dynamic_cast,
//       // since I could keep a separate vector of tabDumpEvents, rather
//       // than filtering for them here.
//       const tabDumpEvent* pTabDumpEvent
// 	= dynamic_cast<const tabDumpEvent*>(pEvent);

//       if(0 != pTabDumpEvent)
// 	{
// 	  pTabDumpEvent->insertTaggedDumpStreamElts(pTaggedDumpStreamsElt);
// 	}
//     }
//   };

  void
  mzrUnit::insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
  {
    // Model elements.
    xmlpp::Element* pModelElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     eltName::model);

    xmlpp::Element* pExplicitSpeciesTagsElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     eltName::explicitSpeciesTags);

    // Give tags for named species.
    std::for_each(speciesByName.begin(),
		  speciesByName.end(),
		  insertExplicitSpeciesTag(pExplicitSpeciesTagsElt));

    // Give all the reactions, using tags to refer to species.
    xmlpp::Element* pTagReactionsElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     eltName::tagReactions);

    // First the reactions that weren't automatically generated.
    std::for_each(userReactions.begin(),
		  userReactions.end(),
		  std::bind2nd(std::mem_fun(&mzrReaction::insertElt),
			       pTagReactionsElt));

    // Now the reactions that were automatically generated.
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
      = utl::dom::mustGetUniqueChild(pModelElt,
				     eltName::time);

    double time = rMolzer.eventQ.getSimTime();
    pTimeElt->set_attribute(eltName::time_secondsAttr,
			    utl::stringify<double>(time));

    // Streams element.
    xmlpp::Element* pStreamsElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     eltName::streams);

    // Emit the species dumped by each of the species streams.
    xmlpp::Element* pTaggedSpeciesStreamsElt
      = pStreamsElt->add_child(eltName::taggedSpeciesStreams);
    std::for_each(speciesStreams.begin(),
		  speciesStreams.end(),
		  std::bind2nd(std::mem_fun
			       (&mzrSpeciesStream::insertDumpedSpeciesTags),
			       pTaggedSpeciesStreamsElt));

    // Emit the species streams dumped in each file.
    xmlpp::Element* pTaggedDumpStreamsElt
      = pStreamsElt->add_child(eltName::taggedDumpStreams);
    std::for_each(dumpStreams.begin(),
		  dumpStreams.end(),
		  std::bind2nd(std::mem_fun
			       (&mzrDumpStream::insertTaggedDumpStreamElts),
			       pTaggedDumpStreamsElt));
  }
}

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

#ifndef CPT_INPUTCAPABILITIES_H
#define CPT_INPUTCAPABILITIES_H

namespace cpt
{
  // This is how a unit indicates what it can parse and parses what it can
  // from a moleculizer-input document.
  class inputCapabilities
  {
    // Names of elements in input files that are handled by a unit.
    //
    // These are connected with the places in the input file where elements may
    // appear that depend on the presence of a particular unit.
    std::set<std::string> modelContentNames;
    std::set<std::string> reactionGenNames;
    std::set<std::string> explicitSpeciesContentNames;
    std::set<std::string> speciesStreamsContentNames;
    std::set<std::string> eventsContentNames;

  public:

    // This allows the parser to just add up the inputCapabilities's of the
    // modules, to form its own inputCapabilities.
    void
    addCap(const inputCapabilities& rCapToAdd)
    {
      modelContentNames.insert
	(rCapToAdd.modelContentNames.begin(),
	 rCapToAdd.modelContentNames.end());

      reactionGenNames.insert
	(rCapToAdd.reactionGenNames.begin(),
	 rCapToAdd.reactionGenNames.end());

      explicitSpeciesContentNames.insert
	(rCapToAdd.explicitSpeciesContentNames.begin(),
	 rCapToAdd.explicitSpeciesContentNames.end());

      speciesStreamsContentNames.insert
	(rCapToAdd.speciesStreamsContentNames.begin(),
	 rCapToAdd.speciesStreamsContentNames.end());

      eventsContentNames.insert
	(rCapToAdd.eventsContentNames.begin(),
	 rCapToAdd.eventsContentNames.end());
    }

    // To build the above sets in unit constructors.
    //
    // Each is connected with a place in the input file schema where elements
    // may appear that depend on the presence of a particular unit.
    void
    addModelContentName(const std::string& rEltName)
    {
      modelContentNames.insert(rEltName);
    }

    void
    addExplicitSpeciesContentName(const std::string& rEltName)
    {
      explicitSpeciesContentNames.insert(rEltName);
    }

    void
    addSpeciesStreamsContentName(const std::string& rEltName)
    {
      speciesStreamsContentNames.insert(rEltName);
    }

    void
    addEventsContentName(const std::string& rEltName)
    {
      eventsContentNames.insert(rEltName);
    }

    void
    addReactionGenName(const std::string& rEltName)
    {
      reactionGenNames.insert(rEltName);
    }

    bool
    handlesModelContentElt(xmlpp::Element* pElt) const
    {
      return (modelContentNames.end()
	      != modelContentNames.find(pElt->get_name()));
    }

  
    // To test if this unit handles a particular reaction generator element.
    bool
    handlesReactionGensContent(xmlpp::Element* pElt) const
    {
      return (reactionGenNames.end()
	      != reactionGenNames.find(pElt->get_name()));
    }

    // To test if this unit handles a particular explictSpecies content element.
    bool
    handlesExplictSpeciesContent(xmlpp::Element* pElt) const
    {
      return (explicitSpeciesContentNames.end()
	      != explicitSpeciesContentNames.find(pElt->get_name()));
    }

    // To test if this unit handles a particular speciesStreams content element.
    bool
    handlesSpeciesStreamsContent(xmlpp::Element* pElt) const
    {
      return (speciesStreamsContentNames.end()
	      != speciesStreamsContentNames.find(pElt->get_name()));
    }

    // To test if this unit handles a particular events content element.
    bool
    handlesEventsContent(xmlpp::Element* pElt) const
    {
      return (eventsContentNames.end()
	      != eventsContentNames.find(pElt->get_name()));
    }
  };
}

#endif // CPT_INPUTCAPABILITIES_H

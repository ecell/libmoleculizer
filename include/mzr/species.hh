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

#ifndef SPECIES_H
#define SPECIES_H

/*! \defgroup speciesGroup Species
  \ingroup mzrGroup
  \brief Species, other state variables, and their families. */

/*! \file species.hh
  \ingroup speciesGroup
  \brief Defines base class for all chemical species. */

#include <algorithm>
#include "domUtils/domUtils.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/stateVar.hh"
#include "mzr/dumpable.hh"

namespace mzr
{
  /*! \ingroup speciesGroup
    \brief Base class for all chemical species. */
  class species :
    public stateVar
  {
  public:
    // For dumping, the number of different species.
    static int speciesCount;

    species(void)
    {
      // Add this species to the tally of all species for species-count output.
      speciesCount++;
    }
  
    virtual
    ~species(void)
    {}

    // For state dump, where species generally are identified by
    // address.
    std::string
    getTag(void) const
    {
      return domUtils::stringify<const species*>(this);
    }

    // For state dump.  This is supposed to return an informative, humanly
    // useful name, but not necessarily unique in this run and certainly
    // not canonical.
    virtual std::string
    getName(void) const
    {
      return getTag();
    }

    /*! \brief Get number of molecules in species. */
    virtual int
    getPop(void) const = 0;

    double
    getConc(mzrUnit& rMzrUnit) const;

    /*! \brief Update population, noting affected reactions.

    In the base class, this doesn't actually do anything to the population
    of the species. The base class method is intended as a "before method" for
    descendant species. */
    virtual void
    update(std::set<reaction*>& rAffectedReactions,
	   int delta);
  };

  // Mixin for molecular weight.
  // 
  // The formula for "arrow" reactions, a la Stochastirator, doesn't use
  // molecular weight.
  class massive
  {
  public:
    virtual ~massive(){}

    virtual double
    getWeight(void) const = 0;
  };

  // Can't go into mzrXcpt.hh, and no other reasonable home for this now.
  class speciesNotMassiveXcpt : public mzrXcpt
  {
    std::string
    mkMsg(xmlpp::Node* pOffendingNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Expected massive species.";
      return msgStream.str();
    }
  public:
    speciesNotMassiveXcpt(xmlpp::Node* pOffendingNode) :
      mzrXcpt(mkMsg(pOffendingNode))
    {}
  };

  massive*
  mustBeMassiveSpecies(xmlpp::Node* pRequestingNode,
		       species* pSpecies)
    throw(speciesNotMassiveXcpt);
}

#endif




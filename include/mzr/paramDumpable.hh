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

#ifndef PARAMDUMPABLE_H
#define PARAMDUMPABLE_H

#include "domUtils/domUtils.hh"
#include "mzr/dumpable.hh"
#include "mzr/paramSpecies.hh"

namespace mzr
{
  // In all of this, the parameter is playing no useful function, and could
  // be replaced with the paramSpecies itself.

  // A query used to decide whether a paramDumpable will dump a given
  // based on its parameter value.
  template<class paramType>
  class paramQuery
  {
  public:
    virtual
    ~paramQuery(void)
    {}

    virtual bool
    operator()(const paramType& rParam) const = 0;
  };

  // paramQuery that always returns true.
  template<class paramType>
  class passAllQuery : public paramQuery<paramType>
  {
  public:
    bool
    operator()(const paramType& rParam) const
    {
      return true;
    }
  };

  // A dumpable that decides whether or not to dump species of which
  // it is notified.
  template<class paramSpeciesType>
  class paramDumpable :
    public speciesDumpable,
    public speciesNotificationTarget<paramSpeciesType>
  {
  public:
    typedef typename paramSpeciesType::param paramType;

  protected:
    static passAllQuery<paramType> passAll;

    const paramQuery<paramType>* pQuery;

    std::vector<paramSpeciesType*> dumpableSpecies;

  public:
    paramDumpable(const std::string& rName,
		  const paramQuery<paramType>* pParamQuery) :
      speciesDumpable(rName),
      pQuery(pParamQuery)
    {}

    // When notified of a new species, use the query to decide
    // whether or not to dump the species.
    void
    notify(paramSpeciesType* pNotifier);

    void
    doDump(std::ostream& rOs) const;

    // For moleculizer-state, where we want to list all the species
    // being dumped by every dumpable.
    void
    insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
      throw(std::exception);
  };

  //////////////////////////////////////////////////////////////////////

  // When notified of a new species, use the query to decide
  // whether or not to dump the species.
  template<class paramSpeciesType>
  void
  paramDumpable<paramSpeciesType>::
  notify(paramSpeciesType* pNewSpecies)
  {
    if((*pQuery)(pNewSpecies->getParam()))
      {
	dumpableSpecies.push_back(pNewSpecies);
      }
  }

  class accumSpeciesPop :
    public std::unary_function<species*, void>
  {
    int& rTot;
  public:
    accumSpeciesPop(int& rTotal) :
      rTot(rTotal)
    {}

    void
    operator()(const species* pSpecies) const
    {
      rTot += pSpecies->getPop();
    }
  };

  template<class paramSpeciesType>
  void
  paramDumpable<paramSpeciesType>::
  doDump(std::ostream& rOs) const
  {
    int tot = 0;
    std::for_each(dumpableSpecies.begin(),
		  dumpableSpecies.end(),
		  accumSpeciesPop(tot));
    rOs << tot;
  }

  class insertSpeciesTag :
    public std::unary_function<species*, void>
  {
    xmlpp::Element* pParentElt;
  public:
    insertSpeciesTag(xmlpp::Element* pParentElement) :
      pParentElt(pParentElement)
    {}

    void
    operator()(const species* pSpecies) const
    {
      xmlpp::Element* pTaggedSpeciesRefElt
	= pParentElt->add_child(eltName::taggedSpeciesRef);

      pTaggedSpeciesRefElt
	->set_attribute(eltName::taggedSpeciesRef_tagAttr,
			pSpecies->getTag());
    }
  };

  template<class paramSpeciesType>
  void
  paramDumpable<paramSpeciesType>::
  insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
    throw(std::exception)
  {
    std::for_each(dumpableSpecies.begin(),
		  dumpableSpecies.end(),
		  insertSpeciesTag(pParentElt));
  }
}

#endif // PARAMDUMPABLE_H

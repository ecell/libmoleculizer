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

#ifndef CPT_MULTIGLOBALSPECIESDUMPABLE_H
#define CPT_MULTIGLOBALSPECIESDUMPABLE_H

#include "utl/dom.hh"
#include "fnd/dumpable.hh"
#include "fnd/sensitive.hh"
#include "fnd/newSpeciesStimulus.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/globalDumpArg.hh"
#include "cpt/stream.hh"

namespace cpt
{
  // These dumpable classes apparently have to be reimplemented for use with
  // globalSpecies, even though the member functions that don't work with
  // globalSpecies don't really need to be instantiated.  g++ basically
  // instantiates these template member functions as part of its way of
  // handling the overall linkage problem with templates.  (For example, the
  // diagnostic does not say where the code was instantiated, as it always
  // does when the template MUST be instantiated.

  // This class is templated on the type of globalSpecies so that it
  // can act as a base class for queryGlobalSpeciesDumpable, reproducing
  // the same relationship as in fnd::querySpeciesDumpable.  This seems like
  // duplication of code, and there might be a better solution.
  template<class globalSpeciesT>
  class multiGlobalSpeciesDumpable :
    public fnd::dumpable<globalDumpArg>,
    public fnd::sensitive<fnd::newSpeciesStimulus<globalSpeciesT> >,
    public cptSpeciesStream
  {
  public:
    typedef globalSpeciesT globalSpeciesType;

  protected:
    std::vector<const globalSpeciesType*> dumpedSpecies;
    
    // For getting compartment names to generate headers.
    //
    // Note that compartmentGraph is available through each globalSpecies, but
    // the set of totalled globalSpecies might be empty.
    const compartmentGraph& rGraph;
    
  public:
    multiGlobalSpeciesDumpable(const std::string& rName,
			       const compartmentGraph& rCompartmentGraph) :
      fnd::dumpable<globalDumpArg>(rName),
      rGraph(rCompartmentGraph)
    {}

    virtual
    void
    respond(const fnd::newSpeciesStimulus<globalSpeciesType>& rStimulus)
    {
      dumpedSpecies.push_back(rStimulus.getSpecies());
    }

    // This has to be overridden (though usefully) so that the
    // virtual function of the same name from fnd::multiSpeciesDumpable
    // doesn't exist. (It has globalSpecies-incompatible code in it.)
    int
    getTotalPop(const globalDumpArg& rDumpArg) const;

    void
    doDump(const globalDumpArg& rDumpArg) const;

    void
    dumpHeader(const globalDumpArg& rDumpArg) const;

    void
    insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
      throw(std::exception);

    void
    insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
      throw(std::exception);
  };
}

#include "cpt/multiGlobalSpeciesDumpableImpl.hh"

#endif // CPT_MULTIGLOBALSPECIESDUMPABLE_H

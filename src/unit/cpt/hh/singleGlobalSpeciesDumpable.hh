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

#ifndef CPT_SINGLEGLOBALSPECIESDUMPABLE_H
#define CPT_SINGLEGLOBALSPECIESDUMPABLE_H

#include "utl/dom.hh"
#include "fnd/dumpable.hh"
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

  class singleGlobalSpeciesDumpable :
    public fnd::dumpable<globalDumpArg>,
    public cptSpeciesStream
  {
    const globalSpecies* pSpecies;
    
  public:
    singleGlobalSpeciesDumpable(const std::string& rName,
				const globalSpecies* pSpeciesToDump) :
      fnd::dumpable<globalDumpArg>(rName),
      pSpecies(pSpeciesToDump)
    {}

    ~singleGlobalSpeciesDumpable(void)
    {}

    const globalSpecies*
    getSpecies(void) const
    {
      return pSpecies;
    }

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

#endif // CPT_SINGLEGLOBALSPECIESDUMPABLE_H

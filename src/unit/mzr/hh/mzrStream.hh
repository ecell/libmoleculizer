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

#ifndef MZRSTREAM_H
#define MZRSTREAM_H

#include <vector>
#include "utl/dom.hh"

namespace mzr
{
  // Interface for state output of species tags for dumpables
  // that dump one or more species populations.
  class mzrSpeciesStream
  {
  public:
    virtual
    ~mzrSpeciesStream(void)
    {}

    // Alas, uses dumpable::getName, so can't actually
    // do this here.
    virtual void
    insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
      throw(std::exception) = 0;

    virtual void
    insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
      throw(std::exception) = 0;
  };

  // Interface for state output for tabDumpEvent.
  class mzrDumpStream
  {
  protected:
    std::vector<const mzrSpeciesStream*> speciesStreams;
    
  public:
    virtual
    ~mzrDumpStream(void)
    {}
    
    void
    addSpeciesStream(const mzrSpeciesStream* pMzrSpeciesStream)
    {
      speciesStreams.push_back(pMzrSpeciesStream);
    }

    virtual void
    insertTaggedDumpStreamElts(xmlpp::Element* pParent) const
      throw(std::exception) = 0;
  };
}

#endif // MZRSTREAM_H

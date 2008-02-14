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

#ifndef PARSEGLOBALSPECIES_H
#define PARSEGLOBALSPECIES_H

#include <vector>
#include "utl/dom.hh"
#include "cpt/compartmentGraph.hh"

namespace cpt
{
  // Routines for parsing out the "core" parts of any global species
  // specification, extracting
  // everything needed to construct a globalSpecies from elements contained in
  // the parent node.
  //
  // These routines are for use as a basis in parsing different
  // kinds of (actual) global species.
  //
  // These should really use const cmlpp::Node*, but there seem to be
  // some const issues in utl::dom.
  std::vector<double>
  parseBoundaryRates(xmlpp::Node* pParentNode,
		     const compartmentGraph& rGraph);
  
  std::vector<int>
  parseCompartmentPops(xmlpp::Node* pParentNode,
		       const compartmentGraph& rGraph);
}

#endif // PARSEGLOBALSPECIES_H

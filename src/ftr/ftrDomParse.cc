//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#include "mol/mzrModMol.hh"
#include "ftr/ftrUnit.hh"
#include "ftr/ftrEltName.hh"
#include "ftr/parseOmniGen.hh"
#include "ftr/parseUniMolGen.hh"

namespace ftr
{
void
ftrUnit::parseDomInput(xmlpp::Element* pRootElement,
xmlpp::Element* pModelElement,
xmlpp::Element* pStreamElt)
throw(std::exception)
{
// This unit only adds a reaction generator for now.
xmlpp::Element* pReactionGensElt
= utl::dom::mustGetUniqueChild(pModelElement,
mzr::eltName::reactionGens);

// Get the omniGen nodes.
xmlpp::Node::NodeList omniGenNodes
= pReactionGensElt->get_children(eltName::omniGen);

// Add omniFam reaction family for each of the generator nodes.
std::for_each(omniGenNodes.begin(),
omniGenNodes.end(),
parseOmniGen(rMzrUnit,
rMolUnit,
rPlexUnit));

// Get the uniMolGen nodes.
xmlpp::Node::NodeList uniMolGenNodes
= pReactionGensElt->get_children(eltName::uniMolGen);

// Add uniMolFam reaction family for each of the generator nodes.
std::for_each(uniMolGenNodes.begin(),
uniMolGenNodes.end(),
parseUniMolGen(rMzrUnit,
rMolUnit,
rPlexUnit));
}
}

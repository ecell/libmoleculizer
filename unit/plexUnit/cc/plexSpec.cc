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

#include <string>
#include <sstream>
#include "mzr/linearHash.hh"
#include "mol/bindingSite.hh"
#include "mol/mol.hh"
#include "mol/molEltName.hh"
#include "plex/plexSpec.hh"
#include "plex/plex.hh"
#include "plex/plexEltName.hh"

namespace plx
{
  xmlpp::Element*
  plexSiteSpec::insertElt(xmlpp::Element* pBindingElt,
			  const plex& rPlex) const throw(std::exception)
  {
    xmlpp::Element* pMolInstanceRefElt
      = pBindingElt->add_child(eltName::molInstanceRef);

    // Cause the mol on which the site occurs to generate the standard
    // fake instance name from the instance index.
    bnd::mol* pMol = rPlex.mols[molNdx()];
    pMolInstanceRefElt->set_attribute(eltName::molInstanceRef_nameAttr,
				      pMol->genInstanceName(molNdx()));

    xmlpp::Element* pBindingSiteRefElt
      = pMolInstanceRefElt->add_child(bnd::eltName::bindingSiteRef);

    // Refer to the binding site by name instead of index.
    const bnd::bindingSite& rBindingSite
      = rPlex.mols[first]->getSiteByNdx(second);
    pBindingSiteRefElt->set_attribute(bnd::eltName::bindingSiteRef_nameAttr,
				      rBindingSite.getName());

    return pMolInstanceRefElt;
  }

  xmlpp::Element*
  plexBinding::insertElt(xmlpp::Element* pPlexElt,
			 const plex& rPlex) const throw(std::exception)
  {
    xmlpp::Element* pBindingElt
      = pPlexElt->add_child(eltName::binding);

    leftSite().insertElt(pBindingElt,
			 rPlex);

    rightSite().insertElt(pBindingElt,
			  rPlex);
    return pBindingElt;
  }
}

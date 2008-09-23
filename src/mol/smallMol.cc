/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008 The Molecular Sciences Institute.
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

#include "mol/smallMol.hh"
#include "mol/molEltName.hh"

namespace bnd
{
  std::vector<mzrBndSite>
  smallMol::
  makeBindingSites(const std::string& rMolName)
  {
      // Heavyweight in order for pedantry.

    const std::string& bindingSiteName( rMolName );

    // Make set of binding site names.
    const std::string& defaultBindingSiteShapeName( rMolName );

      std::set<std::string> bindingSiteShapeNames;
      bindingSiteShapeNames.insert(defaultBindingSiteShapeName);
      
    // Make the single binding site.
    mzrBndSite theSite(bindingSiteName,
		       bindingSiteShapeNames,
		       defaultBindingSiteShapeName);

    // Return vector just containing the one binding site.
    return std::vector<mzrBndSite>(1, theSite);
  }

  std::string
  smallMol::
  genInstanceName(int molInstanceNdx) const
  {
    std::ostringstream oss;
    oss << "small-mol_"
	<< molInstanceNdx;
    return oss.str();
  }

  // Still have to amend state schema to include these elements.
  xmlpp::Element* 
  smallMol::insertElt(xmlpp::Element* pMolsElt) const throw(std::exception)
  {
    // Insert the head element for this smallMol.
    xmlpp::Element* pSmallMolElt
      = pMolsElt->add_child(eltName::smallMol);

    // Add the smallMol's name to the head element.
    pSmallMolElt->set_attribute(eltName::smallMol_nameAttr,
				getName());

    // Insert the weight element.
    xmlpp::Element* pWeightElt
      = pSmallMolElt->add_child(eltName::weight);

    // Add the mol weight in the daltons attribute.
    const cpx::molState* pDefaultState 
      = getDefaultState();
    double molWeight 
      = pDefaultState->getMolWeight();
    pWeightElt->set_attribute(eltName::weight_daltonsAttr,
			      utl::stringify<double>(molWeight));

    return pSmallMolElt;
  }
}

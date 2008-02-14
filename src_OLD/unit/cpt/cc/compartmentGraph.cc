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

#include <numeric>
#include "cpt/cptEltName.hh"
#include "cpt/compartmentGraph.hh"
#include "cpt/unkBoundaryXcpt.hh"
#include "cpt/dupConnectionXcpt.hh"
#include "cpt/dupCompartmentNameXcpt.hh"
#include "cpt/unkCompartmentXcpt.hh"

namespace cpt
{
  void
  compartmentGraph::
  mustAddCompartment(compartment* pNewCompartment,
		     xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    int newCompartmentIndex = compartments.size();
    
    std::pair<cptMapType::iterator, bool> insertResult
      = compartmentsByName.insert
      (cptMapType::value_type(pNewCompartment->getName(),
			      newCompartmentIndex));
    if(insertResult.second)
      {
	compartments.push_back(pNewCompartment);
	pNewCompartment->setIndex(newCompartmentIndex);
      }
    else
      {
	throw dupCompartmentNameXcpt(pNewCompartment->getName(),
				     pRequestingNode);
      }
  }

  bool
  compartmentGraph::
  findCompartmentIndex(const std::string& rCompartmentName,
		       int& rCompartmentIndex) const
  {
    cptMapType::const_iterator iNdx
      = compartmentsByName.find(rCompartmentName);

    if(compartmentsByName.end() == iNdx)
      {
	return false;
      }
    else
      {
	rCompartmentIndex = iNdx->second;
	return true;
      }
  }

  int
  compartmentGraph::
  mustFindCompartmentIndex(const std::string& rCompartmentName,
			   xmlpp::Node* pReferringNode) const
    throw(utl::xcpt)
  {
    int compartmentIndex = -1;
    if(findCompartmentIndex(rCompartmentName,
			    compartmentIndex))
      {
	return compartmentIndex;
      }
    else
      {
	throw unkCompartmentXcpt(rCompartmentName,
				 pReferringNode);
      }
  }

  void
  compartmentGraph::
  mustAddBoundary(boundary* pNewBoundary,
		  xmlpp::Node* pNewBoundaryNode)
    throw(utl::xcpt)
  {
    // Is there already a boundary between the compartments?
    //
    // We can't do the "insert first" trick here, since the map
    // from pairs of components to boundaries
    int boundaryIndex = -1;
    if(findBoundaryIndex(pNewBoundary->getFirstCpt(),
			 pNewBoundary->getSecondCpt(),
			 boundaryIndex))
      {
	throw dupConnectionXcpt(pNewBoundary,
				boundaries[boundaryIndex],
				pNewBoundaryNode);
      }
    else
      {
	boundaryIndex = boundaries.size();

	boundaries.push_back(pNewBoundary);

	// Insert the boundary into the index with compartments
	// in both orders, so search can look in arbitrary order.
	boundariesByCpt.insert
	  (bndMapType::value_type
	   (bndMapType::key_type(pNewBoundary->getFirstCpt(),
				 pNewBoundary->getSecondCpt()),
	    boundaryIndex));

	boundariesByCpt.insert
	  (bndMapType::value_type
	   (bndMapType::key_type(pNewBoundary->getSecondCpt(),
				 pNewBoundary->getFirstCpt()),
	    boundaryIndex));
      }
  }

  bool
  compartmentGraph::
  findBoundaryIndex(compartment* pLeftCompartment,
		    compartment* pRightCompartment,
		    int& rBoundaryIndex) const
  {
    // Boundaries are indexed on their compartments in both orders,
    // so we can search with them in arbitrary order.
    bndMapType::const_iterator iNdx
      = boundariesByCpt.find(bndMapType::key_type(pLeftCompartment,
						  pRightCompartment));
    if(boundariesByCpt.end() == iNdx)
      {
	return false;
      }
    else
      {
	rBoundaryIndex = iNdx->second;
	return true;
      }
  }

  int
  compartmentGraph::
  mustFindBoundaryIndex(compartment* pLeftCompartment,
			compartment* pRightCompartment,
			xmlpp::Node* pReferringNode) const
    throw(utl::xcpt)
  {
    int boundaryIndex = -1;
    if(findBoundaryIndex(pLeftCompartment,
			 pRightCompartment,
			 boundaryIndex))
      {
	return boundaryIndex;
      }
    else
      {
	throw unkBoundaryXcpt(pLeftCompartment,
			      pRightCompartment,
			      pReferringNode);
      }
  }

  boundary*
  compartmentGraph::
  findBoundary(compartment* pLeftCompartment,
	       compartment* pRightCompartment) const
  {
    int boundaryIndex = -1;
    return (findBoundaryIndex(pLeftCompartment,
			      pRightCompartment,
			      boundaryIndex))
      ? boundaries[boundaryIndex]
      : 0;
  }

  double
  compartmentGraph::
  getTotalVolume(void) const
  {
    std::vector<double> compartmentVolumes(compartments.size());

    std::transform(compartments.begin(),
		   compartments.end(),
		   compartmentVolumes.begin(),
		   std::mem_fun(&compartment::getVolume));

    return std::accumulate(compartmentVolumes.begin(),
			   compartmentVolumes.end(),
			   0.0);
  }

  void
  compartmentGraph::
  insertState(xmlpp::Element* pModelElt) const
    throw(std::exception)
  {
    xmlpp::Element* pCompartmentGraphElt
      = pModelElt->add_child(eltName::compartmentGraph);

    xmlpp::Element* pCompartmentsElt
      = pCompartmentGraphElt->add_child(eltName::compartments);

    std::for_each(compartments.begin(),
		  compartments.end(),
		  std::bind2nd(std::mem_fun(&compartment::insertState),
			       pCompartmentsElt));

    xmlpp::Element* pBoundariesElt
      = pCompartmentGraphElt->add_child(eltName::boundaries);

    std::for_each(boundaries.begin(),
		  boundaries.end(),
		  std::bind2nd(std::mem_fun(&boundary::insertState),
			       pBoundariesElt));
  }
}

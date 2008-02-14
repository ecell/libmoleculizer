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

#ifndef COMPARTMENTGRAPH_H
#define COMPARTMENTGRAPH_H

#include "utl/autoVector.hh"
#include "utl/dom.hh"
#include "cpt/compartment.hh"
#include "cpt/boundary.hh"

namespace cpt
{
  class compartmentGraph
  {
    // Map allowing one to look up boundaries between compartments.
    //
    // Each boundary index occurs as the value of two pairs of compartments,
    // the compartments adjacent to the boundary, listed in both orders.
    typedef
    std::map<std::pair<compartment*, compartment*>, int>
    bndMapType;

    bndMapType boundariesByCpt;

  public:
    // For parsing.  This could go into an auxiliary parsing class,
    // similar to parserPlex.
    typedef std::map<std::string, int> cptMapType;
    cptMapType compartmentsByName;

    // The indexing of the compartments in this vector should be paralleled
    // by the indexing of the populations in each species.
    utl::autoVector<compartment> compartments;

    // The boundaries between comparments, which are edges in the compartment
    // graph.
    utl::autoVector<boundary> boundaries;

    // Throws an error if the new compartment's name duplicates
    // any existing compartment name.
    void
    mustAddCompartment(compartment* pNewCompartment,
		       xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // Returns false if no such compartment.
    bool
    findCompartmentIndex(const std::string& rCompartmentName,
			 int& rCompartmentIndex) const;

    // Throws excepton if there is no compartment with
    // the given name.
    int
    mustFindCompartmentIndex(const std::string& rCompartmentName,
			     xmlpp::Node* pReferringNode = 0) const
      throw(utl::xcpt);

    // Returns null if no compartment with the given name.
    compartment*
    findCompartment(const std::string& rCompartmentName) const
    {
      int compartmentIndex = -1;

      return findCompartmentIndex(rCompartmentName,
				  compartmentIndex)
	? compartments[compartmentIndex]
	: 0;
    }

    // Throws an exception if no compartment with the given name.
    compartment*
    mustFindCompartment(const std::string& rCompartmentName,
			xmlpp::Node* pReferringNode = 0) const
      throw(utl::xcpt)
    {
      return
	compartments[mustFindCompartmentIndex(rCompartmentName,
					      pReferringNode)];
    }

    // Adds a boundary component between two compartments (i.e. adds an edge
    // to the compartment graph.)  Throws an exception if the two compartments
    // already have a boundary component between them.
    void
    mustAddBoundary(boundary* pNewBoundary,
		    xmlpp::Node* pNewBoundaryNode = 0)
      throw(utl::xcpt);


    // Returns false if there is no boundary between the given
    // compartments.  Order of the compartments in the calling sequence
    // is arbitrary.
    bool
    findBoundaryIndex(compartment* pLeftCompartment,
		      compartment* pRightCompartment,
		      int& rBoundaryIndex) const;

    // Throws an execption if no boundary between the given
    // compartments.  Order of the compartments is arbitrary.
    int
    mustFindBoundaryIndex(compartment* pLeftCompartment,
			  compartment* pRightCompartment,
			  xmlpp::Node* pReferringNode = 0) const
      throw(utl::xcpt);

    // Returns null if no boundary between the given compartments.
    // Order of the compartments is arbitrary.
    boundary*
    findBoundary(compartment* pLeftCompartment,
		 compartment* pRightCompartment) const;

    // Throws an exception if no boundary between the given compartments.
    // Compartment order is arbitrary.
    boundary*
    mustFindBoundary(compartment* pLeftCompartment,
		     compartment* pRightCompartment,
		     xmlpp::Node* pReferringNode = 0) const
      throw(utl::xcpt)
    {
      return boundaries[mustFindBoundaryIndex(pLeftCompartment,
					      pRightCompartment,
					      pReferringNode)];
    }

    double
    getTotalVolume(void) const;

    void
    insertState(xmlpp::Element* pModelElt) const
      throw(std::exception);
  };
}

#endif // COMPARTMENTGRAPH_H

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

#ifndef PLEXDOMPARSE_H
#define PLEXDOMPARSE_H

#include "domUtils/domUtils.hh"
#include "mol/molState.hh"
#include "plex/prm.hh"
#include "plex/plex.hh"
#include "plex/plexUnit.hh"

namespace plx
{
  // This is also used, e.g. in gpaExchange parsing, to parse the
  // enabling complex.
  plexFamily*
  recognizePlexElt(xmlpp::Element* pPlexElt,
		   parserPlex& rParserPlex,
		   bnd::molUnit& rMolUnit,
		   plexUnit& rPlexUnit);

  // Replaces (default) molParams with those given by modification maps
  // for the corresponding mod-mol instances.
  //
  // This is also used, e.g. in gpaExchange parsing.
  class replaceModMolDefaultState :
    public std::unary_function<const xmlpp::Node*, void>
  {
    std::vector<bnd::molParam>& rParams;
    const plexFamily* pFamily;
    const parserPlex& rPlex;
    bnd::molUnit& rMolUnit;
  
  public:
    replaceModMolDefaultState(std::vector<bnd::molParam>& rMolParams,
			      const plexFamily* pPlexFamily,
			      const parserPlex& rParserPlex,
			      bnd::molUnit& refMolUnit) :
      rParams(rMolParams),
      pFamily(pPlexFamily),
      rPlex(rParserPlex),
      rMolUnit(refMolUnit)
    {}

    void
    operator()(const xmlpp::Node* pModMolInstanceRefNode) const
      throw(mzr::mzrXcpt);
  };

  // Class for general parsing of a plex class (plex with state
  // specifications).  This could always have been used in parsing omni
  // dumpables and plex dumpables. It can now be used in in allosteric-plex
  // and allosteric-omni parsing, as well as in for the new bndKinase reaction
  // generator.
  //
  // It produces the recognized complex, together with a query which, when
  // applied to a species of the recognized complex, tells whether the species
  // satisfies the state specifications.
  class parsePlexClass : public
  std::unary_function<xmlpp::Node*, plexFamily*>
  {
    // For "database access."
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;

    // Output variables.
    parserPlex& rPlex;
    andPlexQueries* pQuery;
  public:
    parsePlexClass(bnd::molUnit& refMolUnit,
		   plexUnit& refPlexUnit,
		   parserPlex& refParsedPlex,
		   andPlexQueries* pParsedQueries) :
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rPlex(refParsedPlex),
      pQuery(pParsedQueries)
    {}

    plexFamily*
    operator()(xmlpp::Node* pParentNode) const;
  };

  // Processes modification map entries into mol state queries, which
  // are added to a query for dumping.
  class addInstanceClausesToQuery :
    public std::unary_function<xmlpp::Node*, void>
  {
    andPlexQueries* pQuery;
    const parserPlex& rPlex;
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    
  public:
    addInstanceClausesToQuery(andPlexQueries* pAndPlexQueries,
			      const parserPlex& rParsedPlex,
			      bnd::molUnit& refMolUnit,
			      plexUnit& refPlexUnit) :
      pQuery(pAndPlexQueries),
      rPlex(rParsedPlex),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}

    void
    operator()(xmlpp::Node* pModMolInstanceRefNode) const
      throw(std::exception);
  };

  // Parses the content of an allosteric-sites element, which conveys a map
  // from binding sites in the complex to binding site shapes.  Each entry in
  // this map comes from parsing a mol-instance-ref element.  The
  // mol-instance-ref element's content conveys a mapping from binding sites
  // on the mol instance to binding site shapes.
  class makeAlloSiteMapEntry : public
  std::unary_function<xmlpp::Node*, std::pair<plexSiteSpec, bnd::siteParam> >
  {
    const parserPlex& rPlex;
  public:
    makeAlloSiteMapEntry(const parserPlex& rParsedPlex) :
      rPlex(rParsedPlex)
    {}

    std::pair<plexSiteSpec, bnd::siteParam>
    operator()(xmlpp::Node* pInstanceRefNode) const
      throw(mzr::mzrXcpt);
  };
}

#endif // PLEXDOMPARSE_H

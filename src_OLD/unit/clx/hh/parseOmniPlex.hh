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

#ifndef CLX_PARSEOMNIPLEX_H
#define CLX_PARSEOMNIPLEX_H

#include "utl/dom.hh"
#include "cml/cmlUnit.hh"
#include "clx/parserPlex.hh"
#include "clx/clxUnit.hh"

namespace clx
{
  // This generates the omniPlex, installs it in its plexFamily, and adds its
  // plexFamily to the clxUnit's list of all omniplex-bearing plexFamilies.
  //
  // Note that the plexFamily is found using "unify," so that if the
  // plexFamily is actually created, because seen for the first time, the
  // plexFamily is not initialized in the usual way
  // (plexFamily::connectToFeatures).  This is because connectToFeatures
  // requires all the omniPlexes to be in place.
  //
  // The way the clxUnit deals with omniPlexes is a two-step process: first,
  // all omniPlexes to be parsed and 'unified into' the database, which must
  // happen before any plexFamily connects to its features.  Second, we sweep
  // through all the omniPlex families again, connecting them to their
  // features.Subsequent "recognitions" can work in the usual way, in which
  // new plexFamilies are initialized (connected to their features)
  // immediately after being created.
  class parseOmniPlex :
    public std::unary_function<xmlpp::Node*, void>
  {
    cpt::cptUnit& rCptUnit;
    cml::cmlUnit& rCmlUnit;
    clxUnit& rClxUnit;
    

  public:
    parseOmniPlex(cpt::cptUnit& refCptUnit,
		  cml::cmlUnit& refCmlUnit,
		  clxUnit& refClxUnit) :
      rCptUnit(refCptUnit),
      rCmlUnit(refCmlUnit),
      rClxUnit(refClxUnit)
    {}

    void
    operator()(xmlpp::Node* pParentNode) const
      throw(utl::xcpt);
  };

  // For use by modules to find an omniPlex parsed for them by the clxUnit.
  //
  // This routine reparses the omniPlex elements below the given parent node,
  // generating a parserPlex.  (This allows module writers to get access to
  // the mol instances, bindings, etc. of the omniPlex for comparison
  // with other things parsed, say, to create a reaction generator.)
  //
  // After findOmni returns, the parserPlex will correspond to the plex below
  // pParentNode, with its mol instances arranged as in the recognizer
  // database; that is, as in the paradigm plex of the omniPlex's plexFamily.
  //
  // With the changes to omniplex implementation, this routine becomes
  // embarrasingly simple, but leaving it as-is for now.
  cptOmniPlex*
  findOmni(xmlpp::Node* pParentNode,
	   cml::cmlUnit& rCmlUnit,
	   clxUnit& rClxUnit,
	   parserPlex& rParsedPlex)
    throw(utl::xcpt);

  // Parses allosteric-omni element.
  class parseAllostericOmni :
    public std::unary_function<xmlpp::Node*, cptOmniPlex*>
  {
    cml::cmlUnit& rCmlUnit;
    clxUnit& rClxUnit;
    
  public:
    parseAllostericOmni(cml::cmlUnit& refCmlUnit,
			clxUnit& refClxUnit) :
      rCmlUnit(refCmlUnit),
      rClxUnit(refClxUnit)
    {}

    cptOmniPlex*
    operator()(xmlpp::Node* pParentNode) const
      throw(utl::xcpt);
  };

  // Parses omni-species-stream element.
  class parseOmniSpeciesStream :
    public std::unary_function<xmlpp::Node*, void>
  {
    cpt::cptUnit& rCptUnit;
    cml::cmlUnit& rCmlUnit;
    clxUnit& rClxUnit;

  public:
    parseOmniSpeciesStream(cpt::cptUnit& refCptUnit,
			   cml::cmlUnit& refCmlUnit,
			   clxUnit& refClxUnit) :
      rCptUnit(refCptUnit),
      rCmlUnit(refCmlUnit),
      rClxUnit(refClxUnit)
    {}

    void
    operator()(xmlpp::Node* pOmniSpeciesStreamNode) const
      throw(utl::xcpt);
  };
}

#endif // CLX_PARSEOMNIPLEX_H

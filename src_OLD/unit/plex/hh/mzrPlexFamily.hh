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

#ifndef PLEXFAMILY_H
#define PLEXFAMILY_H

/*! \defgroup plexStructGroup Structure
  \ingroup plexGroup
  \brief Structural equivalence, structural families of complexes. */

/*! \file mzrPlexFamily.hh
  \ingroup plexStructGroup
  \brief Defines plexFamily, a structural family of species of complexes. */

#include "cpx/plexFamily.hh"
#include "plex/mzrPlex.hh"
#include "plex/mzrPlexSpecies.hh"

namespace plx
{
  class mzrOmniPlex;
  
  /*! \ingroup plexStructGroup

  \brief A structural family of species of complexes.

  The parameter that is used to classify complexes is the vector
  of molParams of the complex.  All other properties of the
  complex are computed from these using the allostery
  function of the structural family. */
  class mzrPlexFamily
    : public cpx::plexFamily<bnd::mzrMol,
			     mzrPlex,
			     mzrPlexSpecies,
			     mzrPlexFamily,
			     mzrOmniPlex>
  {
  public:
    // The arguments other than the paradigm plex are passed on to the base
    // class constructor.  The knownBindings and the set of all omniPlexes
    // are maintained by the plexUnit.
    mzrPlexFamily(const mzrPlex& rParadigm,
		  cpx::knownBindings<bnd::mzrMol, fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> > >& refKnownBindings,
		  std::set<mzrPlexFamily*>& refOmniplexFamilies);

    // Fulfills plexFamily::constructSpecies pure virtual function.
    // This exists so that mzrPlexSpecies can have a reference to
    // the precise mzrPlexFamily class.
    mzrPlexSpecies*
    constructSpecies(const cpx::siteToShapeMap& rSiteParams,
		     const std::vector<cpx::molParam>& rMolParams);

    // Output routine.
    void
    insertSpecies(xmlpp::Element* pExplicitSpeciesElt,
		  double molarFactor) const
      throw(std::exception);
  };
}

#endif


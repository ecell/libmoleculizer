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

#ifndef UNITSMGR_H
#define UNITSMGR_H

#include "mzr/util.hh"
#include "mzr/unit.hh"

#include "mzr/mzrUnit.hh"
#include "mol/molUnit.hh"
#include "dimer/dimerUnit.hh"
#include "plex/plexUnit.hh"
#include "stoch/stochUnit.hh"
#include "gpa/gpaUnit.hh"
#include "nucEx/nucExUnit.hh"
#include "modKinase/modKinaseUnit.hh"
#include "scaffold/scaffoldUnit.hh"
#include "bndKinase/bndKinaseUnit.hh"

namespace mzr
{
  // As a vector, the units should be ordered in output order, so that if each
  // unit arranges its elements correctly (under various heads in the output
  // document) and the units are ordered thus, then all the elements in the
  // output doc should appear in the order demanded by the schema.
  class unitsMgr : public std::vector<unit*>
  {

    // Class for constructing overall input capabilities of moleculizer
    // from input capabilities of each unit.
    class addCapToUnion :
      public std::unary_function<const unit*, void>
    {
      inputCapabilities& rUnionCaps;
    public:
      addCapToUnion(inputCapabilities& rUnionCapabilities) :
	rUnionCaps(rUnionCapabilities)
      {}
    
      void
      operator()(const unit* pUnit) const
      {
	rUnionCaps.addCap(pUnit->inputCap);
      }
    };
      
  public:

    mzr::mzrUnit theMzrUnit;
    bnd::molUnit theMolUnit;
    plx::plexUnit thePlexUnit;
    dimer::dimerUnit theDimerUnit;
    stoch::stochUnit theStochUnit;
//     gpa::gpaUnit theGpaUnit;
//     nucEx::nucExUnit theNucExUnit;
//     kinase::modKinaseUnit theModKinaseUnit;
//     scaf::scaffoldUnit theScaffoldUnit;
    bndKinase::bndKinaseUnit theBndKinaseUnit;

    unitsMgr(moleculizer& rMoleculizer) :
      theMzrUnit(rMoleculizer),
      theMolUnit(rMoleculizer,
		 theMzrUnit),
      thePlexUnit(rMoleculizer,
		  theMzrUnit,
		  theMolUnit),
      theDimerUnit(rMoleculizer,
		   theMzrUnit,
		   theMolUnit,
		   thePlexUnit),
      theStochUnit(rMoleculizer,
		   theMzrUnit),
//       theGpaUnit(rMoleculizer,
// 		 theMzrUnit,
// 		 theMolUnit,
// 		 thePlexUnit,
// 		 theStochUnit),
//       theNucExUnit(rMoleculizer,
// 		   theMzrUnit,
// 		   theMolUnit,
// 		   thePlexUnit,
// 		   theStochUnit),
//       theModKinaseUnit(rMoleculizer,
// 		       theMzrUnit,
// 		       theMolUnit,
// 		       thePlexUnit,
// 		       theStochUnit),
//       theScaffoldUnit(rMoleculizer,
// 		      theMzrUnit,
// 		      theMolUnit,
// 		      thePlexUnit),
      theBndKinaseUnit(rMoleculizer,
		       theMzrUnit,
		       theMolUnit,
		       thePlexUnit)
    {
      // Note that these need to be in output order, which is slightly
      // different from linkage order: mzrUnit "plays cleanup" by parsing
      // after other units are done.
      push_back(&theMolUnit);
      push_back(&theDimerUnit);
      push_back(&thePlexUnit);
      push_back(&theStochUnit);
//       push_back(&theGpaUnit);
//       push_back(&theNucExUnit);
//       push_back(&theModKinaseUnit);
//       push_back(&theScaffoldUnit);
      push_back(&theMzrUnit);
      push_back(&theBndKinaseUnit);
    }

    // Prepares the overall input capabilities of moleculizer from
    // the several input capabilities of the units.
    void
    unionInputCaps(inputCapabilities& rUnion)
    {
      for_each(begin(),
	       end(),
	       addCapToUnion(rUnion));
    }
  };
}

#endif // UNITSMGR_H

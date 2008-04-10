/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008  Nathan Addy
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
//   Nathan Addy, Research Associate     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef __COMPLEXSPECIES_HH
#define __COMPLEXSPECIES_HH

#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <map>

#include "csUtl.hh"
#include "csException.hh"
#include "permutation.hh"

#include "complexOutputState.hh"
#include "partialTokenList.hh"

namespace nmr
{

  // Rewriting ComplexSpecies to operate in terms of Mol*'s.
  class ComplexSpecies
  {

  public:
    typedef std::string Alias;
    typedef int MolNdx;
    typedef std::vector<Mol*> MolList;
    typedef typename Mol::BindingSite BindingSite;
    typedef typename Mol::ModificationList ModificationList;
    typedef int BndNdx;
    typedef std::pair<std::pair<MolNdx, BndNdx>, std::pair<MolNdx, BndNdx> > Binding;
    typedef std::vector<Binding> BindingList;

  public:

    ComplexSpecies();
    ComplexSpecies(const ComplexSpecies& aComplexSpecies);
    ComplexSpecies(const detail::ComplexOutputState& aComplexOutputState)
    {
      // TODO: Write me, as I am used in nmrUnit::getSpeciesFromName.
    }
    ~ComplexSpecies();

    void addMolToComplex(Mol* ptrMol, Alias anAlias);
    void addBindingToComplex(Alias firstMolAlias, BindingSite firstMolSiteAlias, Alias secondMolAlias, BindingSite secondMolSiteAlias);
    int getNumberOfMolsInComplex() const;
    int getNumberOfBindingsInComplex() const;
  
    const MolList& getMolList() const;
    const BindingList& getBindingList() const;
    void applyPermutationToComplex(detail::Permutation& aPermutation);

    // This is a partialNameSentence with everything "stringified".  This plexOutputState is 
    // the object passed to a particular nameAssembler.
    void constructOutputState(detail::ComplexOutputState& rOutputState) const;

  protected:

    typedef std::map<Alias, MolNdx> MolMap;
    typedef MolMap::iterator MolMapIter;

    typedef detail::ComplexOutputState::MolTokenStr MolTokenStr;
    typedef detail::ComplexOutputState::BindingTokenStr BindingTokenStr;
    typedef detail::ComplexOutputState::ModificationTokenStr ModificationTokenStr;

    void constructPartialTokenList(detail::PartialTokenList<molT>& rComplexPartialTokenList) const;
    void sortBinding(Binding& aBinding);

    MolMap theMolAliasToNdxMap;  
    MolList theMols;
    BindingList theBindings;

  };

}

#endif


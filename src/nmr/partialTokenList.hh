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


#ifndef __PARTIALTOKENLIST_HH
#define __PARTIALTOKENLIST_HH

#include <vector>
#include <utility>
#include <iterator>

namespace nmr
{

  namespace detail
  {

    template <class molT>
    struct PartialTokenList
    {
      // PartialTokenLists are designed to compare isomorphic lists -- ie lists produced by indentical
      // complexes.  At least operator< and operator== use this assumption in their code.
      // Do not try to compare non-isomorphic PartialTokenList's.  Results will be undefined, and
      // will likely result in a crash.

      typedef molT Mol;
      typedef int MolNdx; 

      typedef int MolBindingSite;    
      typedef std::pair<MolNdx, MolBindingSite> BindingSite;
      typedef std::pair<BindingSite, BindingSite> Binding;

      typedef std::string ModificationSite;
      typedef std::string ModificationValue;
      typedef std::pair<MolNdx, std::pair<ModificationSite, ModificationValue> > Modification;


      typedef std::vector<Mol> MolList;
      typedef std::vector<Binding> BindingList;
      typedef BindingList::iterator BindingListIter;
      typedef BindingList::const_iterator ConstBindingListIter;
      typedef std::vector<Modification> ModificationList;
      typedef ModificationList::iterator ModificationListIter;
      typedef ModificationList::const_iterator ConstModificationListIter;

      PartialTokenList()
	:
	theMols(),
	theBindings(),
	theModifications(),
	isComplete()
      {
	; // do nothing
      }

      PartialTokenList(const PartialTokenList<molT>& aPns) : theMols(aPns.theMols.begin(), aPns.theMols.end()),
							     theBindings(aPns.theBindings.begin(), aPns.theBindings.end()),
							     theModifications(aPns.theModifications.begin(), aPns.theModifications.end()),
							     isComplete(aPns.isComplete)
      {}

      void clear()
      {
	theMols.clear();
	theBindings.clear();
	theModifications.clear();
      }

    
      MolList theMols;
      BindingList theBindings;
      ModificationList theModifications;  
      bool isComplete;      
      
      bool isSubsetOf(const PartialTokenList<molT>& aPns) const;
      bool isEquivalentTo(const PartialTokenList<molT>& aPns) const;
      bool operator<(const PartialTokenList<molT>& aPns) const; 
      bool operator==(const PartialTokenList<molT>& aPns) const;

    };

    template <class molT>
    bool PartialTokenList<molT>::isSubsetOf(const PartialTokenList<molT>& aPns) const
    {
      if (theBindings.size() > aPns.theBindings.size())
	{
	  return false;
	}
      
      ConstBindingListIter thisBindingListIter;
      ConstBindingListIter aPnsBindingListIter;
      

      for(thisBindingListIter = theBindings.begin(), aPnsBindingListIter = aPns.theBindings.begin();
	  thisBindingListIter != theBindings.end();
	  ++thisBindingListIter, ++aPnsBindingListIter)
	{
	  if (*thisBindingListIter != *aPnsBindingListIter)
	    {
	      return false;
	    }
	}

      // If we make it here, the binding list of this is a subset of the binding list of aPns.
      // Therefore aPns does not have fewer modifications in its ModificationList than this.
      // (If this has an incomplete binding list, it will have 0 modifications in it's list.
      // If it has a complete binding list, then so does aPns, which means that every modification
      // in the complex being studied is in the aPns MolList.
      
      ConstModificationListIter thisModificationListIter;
      ConstModificationListIter aPnsModificationListIter;
      
      for(thisModificationListIter = theModifications.begin(), aPnsModificationListIter = aPns.theModifications.begin();
	  thisModificationListIter != theModifications.end();
	  ++thisModificationListIter, ++aPnsModificationListIter)
	{
	  if (*thisModificationListIter != *aPnsModificationListIter)
	    {
	      return false;
	    }
	}

      return true;
    }

    template <class molT>
    bool PartialTokenList<molT>::operator<(const PartialTokenList<molT>& aPns) const
    {
      // Important! This function is only properly defined when this and aPns come from 
      // the some complex.  Don't try to compare PartialTokenLists from two 
      // non-isomorphic complexes.  It may or may not crash, and it defenitly won't give you 
      // meaningful results.


      // Beacuse we assume that both aPns and this come from trying to canonicalize the same
      // complex, we can assume their MolLists are identical, and therefore we start our comparisons
      // with the binding lists.
      
      ConstBindingListIter thisBindingListIter;
      ConstBindingListIter aPnsBindingListIter;

      for(thisBindingListIter = theBindings.begin(), aPnsBindingListIter = aPns.theBindings.begin();
	  thisBindingListIter != theBindings.end() && aPnsBindingListIter != aPns.theBindings.end();
	  ++thisBindingListIter, ++aPnsBindingListIter)
	{
	  if (*thisBindingListIter < *aPnsBindingListIter)
	    {
	      return true;
	    }
	  if (*aPnsBindingListIter > *thisBindingListIter)
	    {
	      return false;
	    }
	}
      
      // If not both are complete, we cannot say.  Therefore return false, for this 
      // cannot be guarenteed to be less than aPns.
      if (!(this->isComplete && aPns.isComplete))
	{
	  return false;
	}

      // So their binding lists are complete, which means that we can completely compare their
      // modification lists.  Since they both come from the same complex, we know they must 
      // have the same number of modifications.

      return theModifications < aPns.theModifications;
    }
      

    template <class molT>
    bool PartialTokenList<molT>::operator==(const PartialTokenList<molT>& aPns) const
    {  
      if (this->isSubsetOf(aPns))
	{
	  return true;
	}
      else if (aPns.isSubsetOf(*this))
	{
	  return true;
	}
      else 
	{
	  return false;
	}
    }

    template <class molT>
    bool PartialTokenList<molT>::isEquivalentTo(const PartialTokenList<molT>& aPns) const
    {
      if (theMols == aPns.theMols && theBindings == aPns.theBindings && theModifications == aPns.theModifications)
	{
	  return true;
	}
      else
	{
	  return false;
	}
    }



  }

}

#endif

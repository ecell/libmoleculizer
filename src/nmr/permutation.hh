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


#ifndef __PERMUTATION_HH
#define __PERMUTATION_HH

#include <ostream>
#include <vector>

#include "csException.hh"

namespace nmr
{

  namespace detail
  {

    const int UNDEF = -1;
  
    struct Permutation
    {

    public:

      Permutation();

      // This makes a completely undefined permutation on n objects.
      Permutation(int n);
      Permutation(const std::vector<int>& aPermutationVector);
      Permutation(const Permutation& aPermutation);
    
      // This constructor copies aPermutation, and then adds the point perm(pos)=value to it.
      Permutation(const Permutation& aPermutation, int pos, int value);
  
      int getValueAtPosition(int pos) const;
      void setValueAtPosition(int pos, int val);
      void resetValueAtPosition(int pos);

      int getValueAtPositionXcpt(int pos) const;
      void setValueAtPositionXcpt(int pos, int val);
      void resetValueAtPositionXcpt(int pos);

      Permutation of(Permutation& compositionPermutation); // Returns the Permutation (*this)( compositionPermutation(x) )
      Permutation invertPermutation();
  
      bool getIsComplete() const;
      bool getIsIncomplete() const ;
      int getPermutationSize() const;
      int& operator[](const int& n); //in the spirit of std::vector, this has no exception throwing version
      int getLeastValueNotInPermutation() const;
      bool checkPermutationLegality();
      bool operator==(const Permutation& pm);
      bool operator<(const Permutation& pm) const;

    protected:
      std::vector<int> thePermutation;
    };

  }
}

#endif


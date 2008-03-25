#include <exception>
#include <functional>
#include <iterator>
#include <algorithm>
#include <map>

#include "nmr/permutation.hh"
#include "nmr/csUtl.hh"

namespace nmr
{

  namespace detail
  {

    Permutation::Permutation()
      :
      thePermutation()
    {
      ; // do nothing 
    }
    
    Permutation::Permutation(int n) 
      : 
      thePermutation(n, UNDEF) //n values of UNDEF
    {
      ; // do nothing
    }

    Permutation::Permutation(const std::vector<int>& aPermutationVector) 
      : 
      thePermutation(aPermutationVector.begin(), aPermutationVector.end())
    {
      ; // do nothing
    }

    Permutation::Permutation(const Permutation& aPermutation) 
      : 
      thePermutation(aPermutation.thePermutation)
    {
      ; // do nothing
    }

    Permutation::Permutation(const Permutation& aPermutation, int pos, int value) 
      : 
    thePermutation(aPermutation.thePermutation)
    {
      //Add some code to throw and exception if the (pos,value) are bad.
      thePermutation[pos]=value;
    }
  
  

    int Permutation::getValueAtPosition(int pos) const
    {
      return this->thePermutation[pos];
    }


    void Permutation::setValueAtPosition(int pos, int val)
    {
      this->thePermutation[pos]=val;
    }

    void Permutation::resetValueAtPosition(int pos)
    {
      this->thePermutation[pos]=UNDEF;
    }

    int Permutation::getValueAtPositionXcpt(int pos) const
    {
      return this->thePermutation.at(pos);
    }

    void Permutation::setValueAtPositionXcpt(int pos, int value)
    {
      this->thePermutation.at(pos);
      this->thePermutation[pos]=value;
    }

    void Permutation::resetValueAtPositionXcpt(int pos)
    {
      this->thePermutation.at(pos);
      this->thePermutation[pos]=UNDEF;
    }

    Permutation Permutation::of(Permutation& compositionPermutation)
    {
      if(!(this->getPermutationSize()==compositionPermutation.getPermutationSize()))
	{
	  throw CSXcpt("Permutation::of(Permutation& compositionPermutation)", "Permutations not the same size");
	}
      Permutation tmpPerm(this->getPermutationSize());
      for(int i=0;i!=this->getPermutationSize();++i)
	{
	  int intermediateValue=compositionPermutation.getValueAtPosition(i);
	  if(intermediateValue==UNDEF)
	    {
	      tmpPerm[i]=UNDEF;
	    }
	  else
	    {
	      int finalValue=this->getValueAtPosition(intermediateValue);
	      tmpPerm[i]=finalValue;
	    }

	}

      return tmpPerm;
    }

    Permutation Permutation::invertPermutation()
    {
      //creates and returns a new Permutation which is the inverse of this one
      Permutation tmpPerm( thePermutation.size());
      std::vector<int>::iterator idx;
      for(unsigned int i=0;i!=thePermutation.size();++i)
	{
	  //for each number, find out if it is in the permutation
	  //if not, continue
	  //if it is, find the index of the permutation that has it as its value and enter
	  //the appropriate entry in tmpPerm
	  idx=find_if(thePermutation.begin(),
		      thePermutation.end(),
		      std::bind1st(std::equal_to<int>(), i));
	  if(idx!=thePermutation.end())
	    {
	      int dist=idx-thePermutation.begin();
	      //so this permutation sends dist->i
	      //thus make i->Dist part of the new perm
	      tmpPerm[i]=dist;
	    }

	
	}
      return tmpPerm;
    }


    bool Permutation::getIsComplete() const
    {
      std::vector<int>::const_iterator i=find(thePermutation.begin(),
					      thePermutation.end(),
					      UNDEF);
      if (i==thePermutation.end())
	return true;
      else return false;
    }


    bool Permutation::getIsIncomplete() const
    {
      return !getIsComplete();
    }

    
    int Permutation::getPermutationSize() const
    {
      return (int) thePermutation.size();
    }

    int& Permutation::operator[](const int& n)
    {
      return this->thePermutation[n];  
    }


    int Permutation::getLeastValueNotInPermutation() const
    {
      //copy the Permutation to a new vector, ommitting any element where the value is less than 0
      //sort the new vector
      //iterate through the new vector, returning the first position such that value!=position
      std::vector<int> positiveValues;
      copy_if(thePermutation.begin(),
	       thePermutation.end(),
	       back_inserter(positiveValues),
	      std::bind2nd(std::greater_equal<int>(), 0));
      std::sort(positiveValues.begin(), 
		positiveValues.end());
      for(int i=0; i!=(int) positiveValues.size();++i)
	{
	  if (positiveValues[i]!=i)
	    {
	      return i;
	    }
	}
      return positiveValues.size();
    }


    bool Permutation::checkPermutationLegality()
    {
      std::vector<int> sizes(this->thePermutation.size(), 0);
      for(std::vector<int>::iterator i=thePermutation.begin();
	  i!=thePermutation.end();
	  ++i)
	{
	  //if (*i) isn't in -1, 0, 1,..., thePermutation.size()-1, return false
	  if( (*i)<UNDEF || static_cast<int>((thePermutation.size()-1)) < (*i))
	    {
	      return 0;
	    }
	  if((*i)!=UNDEF)
	    {
	      sizes[*i]+=1;
	    }
	}
      for(std::vector<int>::iterator i=sizes.begin();
	  i!=sizes.end();
	  ++i)
	{
	  if ((*i)>1)
	    {
	      return 0;
	    }
	}
      return 1;
    }


    bool Permutation::operator==(const Permutation& pm)
    {
      if ((this->thePermutation)==(pm.thePermutation))
	return true;
      else return false;
    }

    bool Permutation::operator<(const Permutation& pm) const
    {
      if (   (this->getPermutationSize())<(pm.getPermutationSize()))
	return true;

      if( (this->getPermutationSize())>(pm.getPermutationSize()))
	return false;

      std::vector<int>::const_iterator i;
      std::vector<int>::const_iterator j;
      j=pm.thePermutation.begin();
      for(i=thePermutation.begin();
	  i!=thePermutation.end();
	  ++i)
	{
	  if (*i<*j) return true;
	  ++j;
	}
      return false;
    }

  }

}

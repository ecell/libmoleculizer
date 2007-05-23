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

#ifndef UTIL_H
#define UTIL_H

/*! \file util.hh
  \ingroup mzrGroup
  \brief Utility containers. */

#include <map>
#include <vector>
#include <algorithm>
#include <string>

namespace mzr
{
  /*! \ingroup mzrGroup
    \brief Function class for deleting map values.
  
    This ham-handed thing should work with multimaps, too? */
  template<class keyType, class valueType>
  class deleteMapValue
    : public std::unary_function<const std::pair<keyType, valueType>&, void>
  {
  public:
    void operator()(const std::pair<keyType,
		    valueType>& rMapPair) const
    {
      delete rMapPair.second;
    }
  };

  /*! \ingroup mzrGroup
    \brief Function class for deleting vector values. */
  template<class valueType>
  class deleteVectorValue :
    public std::unary_function<valueType, void>
  {
  public:
    void operator()(valueType value) const
    {
      delete value;
    }
  };

  /*! \ingroup mzrGroup
    \brief A string-keyed database of pointers to objects. */
  template<class catObject>
  class catalog
    : public std::map<std::string, catObject*>
  {
  public:

    typedef typename std::map<std::string, catObject*>::const_iterator const_iterator;
    typedef typename std::map<std::string, catObject*>::iterator iterator;
  
    virtual 
    ~catalog(void)
    {}
  
    catObject*
    findEntry(const std::string& rObjectName) const
    {
      const_iterator iNamePtr;
      iNamePtr = this->find(rObjectName);
      return this->end() == iNamePtr
	? 0
	: iNamePtr->second;
    }

    /*!  Returns false if the insert fails because there is already an
      entry with the given name. */
    bool
    addEntry(const std::string& rObjectName,
	     catObject* pObject)
    {
      std::pair<iterator, bool> insertResult
	= insert(std::pair<std::string, catObject*>(rObjectName,
						    pObject));
      return insertResult.second;
    }
  };

  /*! \ingroup mzrGroup
    \brief Adds automatic destruction of the values to catalog. */
  template<class catObject>
  class autoCatalog
    : public catalog<catObject>
  {
    class doDelete
      : public std::unary_function<std::pair<const std::string, catObject*>, void>
    {
    public:
      void operator()(std::pair<const std::string, catObject*>& rNamePtrPair) const
      {
	delete rNamePtrPair.second;
      }
    };

  public:
    ~autoCatalog(void)
    {
      for_each(this->begin(),
	       this->end(),
	       doDelete());
    }
  };

  /*! \ingroup mzrGroup
    \brief A vector of pointers that automatically deletes the pointers. */
  template<class autoObject>
  class autoVector
    : public std::vector<autoObject*>
  {
    class doDelete
      : public std::unary_function<autoObject*, void>
    {
    public:
      void operator()(autoObject* pObject) const
      {
	delete pObject;
      }
    };

  public:
    ~autoVector(void)
    {
      for_each(this->begin(),
	       this->end(),
	       doDelete());
    }

    void
    addEntry(autoObject* pObject)
    {
      push_back(pObject);
    }
  };

  /*! \ingroup mzrGroup
    \brief Create new or overwrite existing value associated with a map key.

    This is the same as

    rTargetMap[key] = value,

    except that the indexing form requires a default constructor for
    values, since it could be used as an lvalue:

    foo = rTargetMap[key]. */
  template<class mapClass>
  void
  forceInsert(mapClass& rTargetMap,
	      const typename mapClass::value_type& rKeyValuePair)
  {
    // Attempt to insert the key/value pair.  This will not succeed if
    // the key is already associated to a value by the map.
    std::pair<typename mapClass::iterator, bool> insertResult
      = rTargetMap.insert(rKeyValuePair);

    bool insertSucceeded = insertResult.second;

    // If the key was already mapped, reset the value associated to it.
    if(! insertSucceeded)
      {
	typename mapClass::iterator iEntry = insertResult.first;
	iEntry->second = rKeyValuePair.second;
      }
  }

  /*! \ingroup mzrGroup
    \brief Auxiliary class to reset multiple map values. */
  template<class mapClass>
  class forceInsertOne :
    public std::unary_function<typename mapClass::value_type, void>
  {
    mapClass& rTarget;
  public:
    forceInsertOne(mapClass& rTargetMap) :
      rTarget(rTargetMap)
    {}
    void operator()(const typename mapClass::value_type& rEntry) const
    {
      forceInsert(rTarget,
		  rEntry);
    }
  };

  /*! \ingroup mzrGroup
    \brief Create new or overwrite existing values associated to many map keys.
  */
  template<class mapClass>
  void
  forceInsert(mapClass& rTargetMap,
	      typename mapClass::const_iterator startIter,
	      typename mapClass::const_iterator stopIter)
  {
    for_each(startIter,
	     stopIter,
	     forceInsertOne<mapClass>(rTargetMap));
  }

  /*! \ingroup mzrGroup
    \brief Auxiliary function class for funcInsert.
  */
  template<class mapClass, class functionClass>
  class funcInsertOne :
    public std::unary_function<typename mapClass::value_type, void>
  {
    mapClass& rTarget;
    const functionClass keyFunction;
  public:
    funcInsertOne(mapClass& rTargetMap,
		  const functionClass& rKeyFunction) :
      rTarget(rTargetMap),
      keyFunction(rKeyFunction)
    {}
    void operator()(const typename mapClass::value_type& rKeyValuePair) const
    {
      forceInsert(rTarget,
		  typename mapClass::value_type(keyFunction(rKeyValuePair.first),
						rKeyValuePair.second));
    }
  };

  /*! \ingroup mzrGroup
    \brief ForceInsert with remapping of keys.

    Function class that applies a given function to the key
    before force-inserting into a given map.  This is used
    to apply a plexMap to the sitesSpecs in (siteSpec, siteParam)
    pairs, where the plexMap is an injection of an omniplex into
    a plex. */
  template<class mapClass, class functionClass>
  void
  funcInsert(mapClass& rTargetMap,
	     const functionClass& rKeyFunction,
	     typename mapClass::const_iterator startIter,
	     typename mapClass::const_iterator stopIter)
  {
    for_each(startIter,
	     stopIter,
	     funcInsertOne<mapClass, functionClass>(rTargetMap,
						    rKeyFunction));
  }

  // For "roll your own" conversion of double to scientific notation, so that
  // it can be formatted in mathML.  The calling sequence is modeled on the
  // library routind frexp and its friends.
  // 
  // The frexp routine does this exaclty right in base 2, and this function is
  // a wrapper for it.  The scientific notation generated by this is not
  // precisely correct, in that the fractional part may not of exactly the
  // right magnitude, but it should be numerically correct anyway.  (The
  // fractional part is derived from the base 2 fractional part without making
  // a detailed correction of its magnitude.)
  double
  frexp10(double num,
	  int& rExponent);
}

#endif

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

#ifndef OBJECTCACHE_H
#define OBJECTCACHE_H

#include <map>
#include <functional>
#include <algorithm>

namespace utl
{
  template<class keyType, class objectType>
  class objectCache :
    public std::map<keyType, objectType*>
  {
  public:
    virtual
    ~objectCache(void)
    {}

    virtual objectType*
    makeMember(const keyType& rKey) = 0;

    // This assumes that the copy constructor for keys is lightweight.
    objectType*
    getMember(const keyType& rKey)
    {
      typename std::pair<typename objectCache::iterator, bool> insertResult
	= insert(typename objectCache::value_type(rKey,
			    (objectType*) 0));

      if(insertResult.second)
	{
	  insertResult.first->second = makeMember(rKey);
	}

      return insertResult.first->second;
    }
  };

  template<class keyType, class objectType>
  class autoCache :
    public objectCache<keyType, objectType>
  {
    class doDelete :
      public std::unary_function<typename autoCache::value_type, void>
    {
    public:
      void
      operator()(const
		 typename doDelete::argument_type&
		 rKeyObjectPtrPair) const
      {
	delete rKeyObjectPtrPair.second;
      }
    };
    
  public:
    virtual
    ~autoCache(void)
    {
      std::for_each(this->begin(),
		    this->end(),
		    doDelete());
    }
  };
}

#endif // OBJECTCACHE_H

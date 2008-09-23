/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 The Molecular Sciences Institute
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
/////////////////////////////////////////////////////////////////////////////


#ifndef ENTITYREGISTRY_HH
#define ENTITYREGISTRY_HH

#include "fnd/macros.hh"

#include <string>
#include <boost/call_traits.hpp>
using boost;

namespace fnd
{

  // This is an abstract class that will describe the interfaces for the two 
  // ReactionNetworkDescription 
  template <typename ENTITY_TYPE, typename ENTITY_KEY>
  class EntityRegistry
  {
  public:
    DECLARE_TYPE( ENTITY_KEY, EntityKey );
    DECLARE_TYPE( ENTITY_TYPE, Entity );
    DECLARE_TYPE( MM_TYPE, MemoryManagementStrategy);

    // 
    virtual ~EntityRegistry()
    {}
    
    // These two should always return successfully, returning true if the functions
    // had to do something to complete them, false otherwise.
    virtual bool RegisterEntity(EntityPtr pEntityType) = 0;
    virtual bool DeregisterEntiy(EntityPtr pEntityType) = 0;
    virtual call_traits<Entity> 

  protected:

    std::map<EntityKey, Entity> theEntityMap;
  };

}

#endif

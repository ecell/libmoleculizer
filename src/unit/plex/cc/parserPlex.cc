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

#include "mzr/moleculizer.hh"
#include "plex/parserPlex.hh"

namespace plx
{
  // This doesn't check for duplicate mol instance names.
  int
  parserPlex::
  addMolByName(const std::string& rName,
	       bnd::mzrMol* pMol)
  {
    int ndx = (int) mols.size();
    mols.push_back(pMol);
    nameToMolNdx[rName] = ndx;
    return ndx;
  }
  
  int
  parserPlex::
  getMolNdxByName(const std::string& rName) const
  {
    std::map<std::string, int>::const_iterator iNameNdxPair
      = nameToMolNdx.find(rName);
      
    return (nameToMolNdx.end() == iNameNdxPair) 
      ? -1
      : iNameNdxPair->second;
  }

  int
  parserPlex::
  mustGetMolNdxByName(xmlpp::Node* pRequestingNode,
		      const std::string& rInstanceName) const
    throw(unkMolInstXcpt)
  {
    int ndx = getMolNdxByName(rInstanceName);
    if(ndx < 0) throw unkMolInstXcpt(rInstanceName,
				     pRequestingNode);
    return ndx;
  }

  bnd::mzrMol*
  parserPlex::getMolByName(const std::string& rName) const
  {
    int molNdx = getMolNdxByName(rName);
    return 0 <= molNdx
      ? mols[molNdx]
      : 0;
  }

  bnd::mzrMol*
  parserPlex::mustGetMolByName(xmlpp::Node* pRequestingNode,
			       const std::string& rInstanceName) const
    throw(unkMolInstXcpt)
  {
    return mols[mustGetMolNdxByName(pRequestingNode,
				    rInstanceName)];
  }
}

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

#include "mzr/linearHash.hh"
#include "plex/plex.hh"
#include "plex/plexMap.hh"
#include "plex/plexEltName.hh"

namespace plx
{
  // This allows constructing the map from sites to bindings on the
  // stack. It is needed when finding the free binding sites.
  //
  // Note that each site is on at most one binding, so we can use a map
  // instead of a multimap.
  void
  plex::makeSiteToBindings(std::map<plexSiteSpec, int>& rSiteToBindings) const
  {
    int bindingNdx = bindings.size();
    while(0 < bindingNdx--)
      {
	const plexBinding& rBinding = bindings[bindingNdx];
	rSiteToBindings.insert(std::make_pair(rBinding.leftSite(), bindingNdx));
	rSiteToBindings.insert(std::make_pair(rBinding.rightSite(), bindingNdx));
      }
  }

  void
  plex::makeFreeSiteVector(std::vector<plexSiteSpec>& rFreeSiteVector) const
  {
    // Make the map from sites to the bindings they're involved in.
    std::map<plexSiteSpec, int> siteToBindings;
    makeSiteToBindings(siteToBindings);

    // Enumerate the sites in this plex.  If a site doesn't appear
    // in the sites to bindings map, it's a free site, so insert it
    // into the vector of free sites.
    int molNdx = mols.size();
    while(0 < molNdx--)
      {
	bnd::mol* pMol = mols[molNdx];
	
	int siteNdx = pMol->getSiteCount();
	while(0 < siteNdx--)
	  {
	    plexSiteSpec theSpec(molNdx, siteNdx);
	    if(siteToBindings.find(theSpec) == siteToBindings.end())
	      rFreeSiteVector.push_back(theSpec);
	  }
      }
  }

  // This allows constructing the map from mols to bindings on the
  // stack. It is needed when finding the connected components of a plex.
  //
  // Note that, unlike sites, mols can be on as many bindings as
  // there are sites on the mol, so we have to use a mulitmap.
  void
  plex::makeMolToBindings(std::multimap<int, int>& rMolToBindings) const
  {
    int bindingNdx = bindings.size();
    while(0 < bindingNdx--)
      {
	const plexBinding& rBinding = bindings[bindingNdx];
	rMolToBindings.insert(std::make_pair(rBinding.leftSite().molNdx(),
					     bindingNdx));
	rMolToBindings.insert(std::make_pair(rBinding.rightSite().molNdx(),
					     bindingNdx));
      }
  }

  void
  plex::pushConnectedBindings(int molNdx,
			      const std::multimap<int, int>& rMolToBindings,
			      plexIsoPair& rIso,
			      plex& component) const
  {
    // Which bindings are on mol at molNdx?
    std::pair<std::multimap<int, int>::const_iterator,
      std::multimap<int, int>::const_iterator> rangeIterators
      = rMolToBindings.equal_range(molNdx);

    // Install all of the connected bindings that aren't already
    // in the connected component.
    for(std::multimap<int, int>::const_iterator iMolBindingNdx = rangeIterators.first;
	rangeIterators.second != iMolBindingNdx;
	iMolBindingNdx++)
      {
	// Has this binding already been done?
	int bindingNdx = iMolBindingNdx->second;
	if(rIso.forward.bindingUnmapped(bindingNdx))
	  {
	    const plexBinding& rBinding = bindings[bindingNdx];

	    // Make sure the mols at the ends of the binding are done.
	    // Rather than try to figure out which is the mol opposite
	    // to molNdx, just checking both of them.
	    int leftMolNdx = rBinding.leftSite().molNdx();
	    int leftTargetNdx = pushConnectedMol(leftMolNdx,
						 rIso,
						 component);
	    int rightMolNdx = rBinding.rightSite().molNdx();
	    int rightTargetNdx = pushConnectedMol(rightMolNdx,
						  rIso,
						  component);

	    // Now put in the binding, with its ends remapped.
	    plexSiteSpec leftSiteSpec(leftTargetNdx,
				      rBinding.leftSite().siteNdx());
	    plexSiteSpec rightSiteSpec(rightTargetNdx,
				       rBinding.rightSite().siteNdx());
	    plexBinding targetBinding(leftSiteSpec,
				      rightSiteSpec);
	    int bindingTargetNdx = component.bindings.size();
	    component.bindings.push_back(targetBinding);
	    rIso.forward.bindingMap[bindingNdx] = bindingTargetNdx;
	    rIso.backward.bindingMap[bindingTargetNdx] = bindingNdx;
	  }
      }
  }

  int
  plex::pushConnectedMol(int molNdx,
			 plexIsoPair& rIso,
			 plex& component) const
  {
    int targetNdx;
    if(rIso.forward.molUnmapped(molNdx))
      {
	targetNdx = (int) component.mols.size();
	component.mols.push_back(mols[molNdx]);
	rIso.forward.molMap[molNdx] = targetNdx;
	rIso.backward.molMap[targetNdx] = molNdx;
      }
    else
      {
	targetNdx = rIso.forward.molMap[molNdx];
      }
    return targetNdx;
  }

  // Loads the connected component of mols[molNdx] in this plex into
  // the plex "component," which should be empty when the function
  // is called.
  //
  // The pointer to a partial isomorphism allows tracking in cases where
  // it's needed.  The partial isomorphism should be as large enough
  // to accomodate the entire complex, which might be connected, and it
  // should be empty.  (see "makeConnectedComponent").
  void
  plex::makeTrackedComponent(int molNdx,
			     plex& component,
			     plexIsoPair& rIso) const
  {
    std::multimap<int, int> molToBindings;
    makeMolToBindings(molToBindings);

    // Push the bindings on this mol to the component (along with
    // their mols).
    pushConnectedMol(molNdx,
		     rIso,
		     component);

    // Propagate over the bindings by "data recursion" on the mols
    // of the connected component.
    for(int compMolNdx = 0;
	compMolNdx < (int) component.mols.size();
	compMolNdx++)
      {
	int srcMolNdx = rIso.backward.molMap[compMolNdx];
	pushConnectedBindings(srcMolNdx,
			      molToBindings,
			      rIso,
			      component);
      }
  }

  void
  plex::makeConnectedComponent(int molNdx,
			       plex& component) const
  {
    plexIsoPair iso(mols.size(),
		    bindings.size());

    makeTrackedComponent(molNdx,
			 component,
			 iso);
  }

  bool
  plex::operator<(const plex& rRightPlex) const
  {
    if(mols < rRightPlex.mols) return true;
    if(rRightPlex.mols < mols) return false;
    return bindings < rRightPlex.bindings;
  }

  class plex::insertBindingElt :
    public std::unary_function<plexBinding, void>
  {
    xmlpp::Element* pPlexElt;
    const plex& rPlx;
  public:
    insertBindingElt(xmlpp::Element* pPlexElement,
		     const plex& rPlex) :
      pPlexElt(pPlexElement),
      rPlx(rPlex)
    {}

    void
    operator()(const plexBinding& rBinding) const throw(std::exception)
    {
      rBinding.insertElt(pPlexElt,
			 rPlx);
    }
  };
			    
  xmlpp::Element*
  plex::insertElt(xmlpp::Element* pParentElt) const throw(std::exception)
  {
    xmlpp::Element* pPlexElt
      = pParentElt->add_child(eltName::plex);
  
    // Insert elements for each mol instance.  I'm using the index
    // to generate instance names.  Users will be disappointed that
    // their own instance names have been forgotten.
    for(int molNdx = 0;
	molNdx < (int) mols.size();
	++molNdx)
      {
	xmlpp::Element* pMolInstanceElt
	  = pPlexElt->add_child(eltName::molInstance);

	bnd::mol* pMol = mols[molNdx];

	// Stringify the instance index as an instance name.
	pMolInstanceElt->set_attribute(eltName::molInstance_nameAttr,
				       pMol->genInstanceName(molNdx));

	xmlpp::Element* pMolRefElt
	  = pMolInstanceElt->add_child(eltName::molRef);

	pMolRefElt->set_attribute(eltName::molRef_nameAttr,
				  mols[molNdx]->getName());
      }

    // Insert elements for each binding.
    std::for_each(bindings.begin(),
		  bindings.end(),
		  insertBindingElt(pPlexElt,
				   *this));

    return pPlexElt;
  }
}

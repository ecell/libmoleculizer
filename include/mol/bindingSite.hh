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

#ifndef BINDINGSITE_H
#define BINDINGSITE_H

#include <vector>
#include <string>
#include <set>
#include "mzr/feature.hh"
#include "mol/molXcpt.hh"
#include "plex/plexSpecies.hh"
#include "plex/prm.hh"

class mol;

namespace bnd
{
  // Dimerization generators connect to a bindingSite as a feature in
  // order to be notified when a new plexSpecies is created in which the
  // binding site is unbound, and therefore available for dimerization.
  //
  // Of course, other reaction generators interested in free binding
  // sites could also listen on this channel.
  class bindingSite :
    public mzr::feature<plx::plexSpecies, plx::plexSiteSpec>
  {
    std::string name;

    // Note that the values of this map are actual site shapes, not
    // pointers.  Must not use an autoCatalog, for example, because
    // autoCatalog doesn't have copy semantics: the copied autoCatalog
    // will delete the pointers when it is deleted, so that the copy
    // will contain invalid pointers.  Moreover, according to
    // the nouvelle regime, site shapes belong to their sites and should
    // go away when their sites go away.
    //
    // The assumption that we can reasonably copy site shapes is quite
    // restrictive on future enhancements (e.g. shape-theoretic) of the
    // site shape notion.
    std::map<std::string, siteShape> shapesByName;

    // Pointer to the usual shape of this site.  I'm going to try naming
    // this in some conventional and brief way.  I don't think that I
    // can maintain the molname:sitename convention, so there will be
    // minor breakage in scripts.
    //
    // An alternative would be to make all binding site shape names
    // explicit, so that the name of the default shape of each site
    // would be given in the mol definition.
    //
    // Note that having this as a pointer makes copying more
    // difficult, the copy's default parameter must be set
    // using a name-based lookup.  An alternative would be to
    // intern the site shapes in a vector, use an index here,
    // and store indices in the name-lookup map.
    siteParam defaultParam;

  public:
    bindingSite(const std::string& rSiteName,
		const std::map<std::string, siteShape>& rNameToShape,
		const std::string& rDefaultShapeName) :
      name(rSiteName),
      shapesByName(rNameToShape)
    {
      setDefaultParamByName(rDefaultShapeName);
    }

    bindingSite(xmlpp::Node* pBindingSiteNode) throw(std::exception);

    class constructor :
      public std::unary_function<xmlpp::Node*, bindingSite>
    {
    public:
      bindingSite
      operator()(xmlpp::Node* pBindingSiteNode) const throw(std::exception)
      {
	return bindingSite(pBindingSiteNode);
      }
    };

    // I need a copy constructor, since the pointer to the default
    // parameter has to be adjusted.  This wasn't the case before when
    // site shapes were permanently interned somewhere else.
    //
    // This is slow and sucky, but binding sites are basically copied
    // only in the mol constructor, when mols are copied,etc.
    bindingSite(const bindingSite& rOriginal) :
      name(rOriginal.name),
      shapesByName(rOriginal.shapesByName)
    {
      setDefaultParamByName(rOriginal.getDefaultParam()->name);
    }

    const std::string&
    getName(void) const
    {
      return name;
    }

    siteParam
    getDefaultParam(void) const
    {
      return defaultParam;
    }

    // This routine fails quietly (because there is no site shape with
    // the given name) by returning null, since it is mainly expected to
    // be invoked in user interaction, and the user interaction routine
    // will need to report the error in a user-intelligible way.
    siteParam
    setDefaultParamByName(const std::string& rShapeName)
    {
      defaultParam = getParam(rShapeName);

      return defaultParam;
    }

    // Exposed to allow traversal of all site shapes, so that default
    // binding kinetics can be created.
    const std::map<std::string, siteShape>&
    getSiteShapeMap(void) const
    {
      return shapesByName;
    }

    // This routine fails quietly (because there is no site shape with
    // the given name) by returning null, since it is mainly expected to
    // be invoked in user interaction, and the user interaction routine
    // will need to report the error in a user-intelligible way.
    siteParam
    getParam(const std::string& rShapeName) const
    {
      std::map<std::string, siteShape>::const_iterator iShape
	= shapesByName.find(rShapeName);

      if(shapesByName.end() == iShape) return 0;
      else return &(iShape->second);
    }

    // Does the same as getParam, i.e. returns the site shape
    // with the given name, or throws an exception
    // if there is no site shape with the given name.
    siteParam
    mustGetParam(xmlpp::Node* pRequestingNode,
		 const mol* pMol,
		 const std::string& rShapeName) const
      throw(unknownSiteShapeXcpt);

    // Returns pointer to interned form of the provided shape.
    // 0 is returned if there is already a shape with this name.
    // This routine fails quietly for the above noted reasons.
    siteParam
    addShape(const std::string& rShapeName,
	     const siteShape& rShape)
    {
      std::pair<std::map<std::string, siteShape>::iterator, bool> insertResult
	= shapesByName.insert(std::make_pair(rShapeName,
					     rShape));
    

      if(insertResult.second) return &(insertResult.first->second);
      else return 0;
    }

    // Returns pointer to interned shape, which is created if it
    // does not exist.  This depends on site shapes being "black
    // boxes" with names.
    siteParam
    internShape(const std::string& rShapeName)
    {
      siteParam param = getParam(rShapeName);
      if(0 == param) return addShape(rShapeName,
				     siteShape(rShapeName));
      else return param;
    }

    xmlpp::Element*
    insertElt(xmlpp::Element* pMolElt) const throw(std::exception);
  };
}

#endif // BINDINGSITE_H

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

#include "mol/molEltName.hh"
#include "mol/modMol.hh"
#include "mol/molXcpt.hh"
#include "plex/plexEltName.hh"

namespace bnd
{
  std::string
  modMol::genInstanceName(int molInstanceNdx) const
  {
    std::ostringstream oss;
    oss << "mod-mol_"
	<< molInstanceNdx;
    return oss.str();
  }

  xmlpp::Element*
  modMol::insertInstanceState(xmlpp::Element* pInstanceStatesElt,
			      int molInstanceNdx,
			      molParam param) const
  {
    // Insert mod-mol-ref/mod-map elements as the description of the
    // modMolState pointed to by param.
    //
    // To save space and be consistent with the input document format, we'd like
    // to omit modification sites that have their default modification (ususally
    // something like "none").

    xmlpp::Element* pModMolInstanceRefElt
      = pInstanceStatesElt->add_child(plx::eltName::modMolInstanceRef);

    pModMolInstanceRefElt
      ->set_attribute(plx::eltName::modMolInstanceRef_nameAttr,
		      genInstanceName(molInstanceNdx));

    // Insert the modification-map element.
    xmlpp::Element* pModMapElt = pModMolInstanceRefElt
      ->add_child(eltName::modMap);

    // Get the state pointed to by param.
    const modMolState& rState = externState(param);
  
    // Get the default state for comparison with the specified state.
    const modMolState& rDefaultState = *(getDefaultState());

    // Insert a mod-site-ref/mod-ref pair for each modification site that
    // is not in its default state.
    //
    // This seems unnecessarily complicated, involving indexes of modifications,
    // or alternatively, marching in parallel through the two vectors of
    // modification pointers and the vector of modification site names;
    // like a ternary version of std::for_each.
    for(int modSiteNdx = 0;
	modSiteNdx < modSiteCount();
	++modSiteNdx)
      {
	// Get pointers to the (interned) actual and default modifications
	// of modification site at modNdx.
	const modification* pActualModification = rState[modSiteNdx];
	const modification* pDefaultModification = rDefaultState[modSiteNdx];

	if(pActualModification != pDefaultModification)
	  {
	    // Insert mod-site-ref element.
	    xmlpp::Element* pModSiteRefElt
	      = pModMapElt->add_child(eltName::modSiteRef);

	    // Add the modification site name attribute.
	    pModSiteRefElt->set_attribute(eltName::modSiteRef_nameAttr,
					  modSiteNames[modSiteNdx]);

	    // Add the mod-ref element, giving the name of the modification.
	    xmlpp::Element* pModRefElt
	      = pModSiteRefElt->add_child(eltName::modRef);

	    pModRefElt->set_attribute(eltName::modRef_nameAttr,
				      pActualModification->name);
	  }
      }
    return pModMolInstanceRefElt;
  }

  class insertSite :
    public std::unary_function<bindingSite, void>
  {
    xmlpp::Element* pMolElt;
  public:
    insertSite(xmlpp::Element* pMolElement) :
      pMolElt(pMolElement)
    {}

    void
    operator()(const bindingSite& rSite) const throw(std::exception)
    {
      rSite.insertElt(pMolElt);
    }
  };

  class insertNonDefault : public
  std::unary_function<std::pair<const modMolState, std::vector<siteParam> >,
			   void>
  {
    xmlpp::Element* pModMolElt;
    const modMol* pMol;
    const modMolState* pDfltState;
    const std::vector<siteParam>& rDefaultShapes;
  public:
    insertNonDefault(xmlpp::Element* pModMolElement,
		     const modMol* pModMol,
		     const modMolState* pDefaultState,
		     const std::vector<siteParam>& rDefaultSiteParams) :
      pModMolElt(pModMolElement),
      pMol(pModMol),
      pDfltState(pDefaultState),
      rDefaultShapes(rDefaultSiteParams)
    {}

    void
    operator()(const argument_type& rAlloMapEntry) const throw(std::exception)
    {
      const std::vector<siteParam>& rEntryShapes
	= rAlloMapEntry.second;

      // We don't want to list states that are not truly allosteric
      // in the sense of having some non-default site shape.
      if(rEntryShapes != rDefaultShapes)
	{
	  // Insert the allosteric-state element.
	  //
	  // Is this element necessary or useful?  It is a wrapper for the
	  // different kinds of state descriptions that different kinds of mols
	  // will have.
	  xmlpp::Element* pAllostericStateElt
	    = pModMolElt->add_child(eltName::allostericState);

	  // Insert the mod-map element, whose content describes the
	  // modification state.
	  xmlpp::Element* pModMapElt
	    = pAllostericStateElt->add_child(eltName::modMap);

	  // For each modification site that is not in its default modification
	  // state, emit a modSiteRef/modRef pair, giving the actual
	  // modification.
	  for(int modSiteNdx = 0;
	      modSiteNdx < pMol->modSiteCount();
	      ++modSiteNdx)
	    {
	      const modMolState* pEntryState = &(rAlloMapEntry.first);

	      const modification* pActualMod = (*pEntryState)[modSiteNdx];
	      const modification* pDefaultMod = (*pDfltState)[modSiteNdx];
	    
	      // Is the modification not the default modification for this site?
	      if(pActualMod != pDefaultMod)
		{
		  xmlpp::Element* pModSiteRefElt
		    = pModMapElt->add_child(eltName::modSiteRef);

		  pModSiteRefElt
		    ->set_attribute(eltName::modSiteRef_nameAttr,
				    pMol->modSiteNames[modSiteNdx]);

		  xmlpp::Element* pModRefElt
		    = pModSiteRefElt->add_child(eltName::modRef);

		  pModRefElt->set_attribute(eltName::modRef_nameAttr,
					    pActualMod->name);
		}
	    }

	  // Insert the site-shape-map element, which gives the
	  // (non-default)shapes of the binding sites when in this allosteric
	  // state.
	  xmlpp::Element* pSiteShapeMapElt
	    = pAllostericStateElt->add_child(eltName::siteShapeMap);

	  // For each binding site that is not in its default shape,
	  // emit a bindingSiteRef/siteShapeRef pair, giving the actual
	  // shape of the binding site.
	  for(int siteNdx = 0;
	      siteNdx < pMol->getSiteCount();
	      ++siteNdx)
	    {

	      const siteShape* pDefaultSiteShape = rDefaultShapes[siteNdx];
	      const siteShape* pActualSiteShape = rEntryShapes[siteNdx];

	      if(pDefaultSiteShape != pActualSiteShape)
		{
		  xmlpp::Element* pBindingSiteRefElt
		    = pSiteShapeMapElt->add_child(eltName::bindingSiteRef);

		  pBindingSiteRefElt
		    ->set_attribute(eltName::bindingSiteRef_nameAttr,
				    pMol->getSiteByNdx(siteNdx).getName());

		  xmlpp::Element* pSiteShapeRefElt
		    = pBindingSiteRefElt->add_child(eltName::siteShapeRef);

		  pSiteShapeRefElt
		    ->set_attribute(eltName::siteShapeRef_nameAttr,
				    pActualSiteShape->name);
		}
	    }
	}
    }
  };

  xmlpp::Element*
  modMol::insertElt(xmlpp::Element* pMolsElt) const throw(std::exception)
  {
    // We need the default state to compute the molecular weight and to
    // get the list of default modifications.
    const modMolState* pDefaultState = getDefaultState();
  
    // Insert the head element for this modMol.
    xmlpp::Element* pModMolElt
      = pMolsElt->add_child(eltName::modMol);

    pModMolElt->set_attribute(eltName::modMol_nameAttr,
			      getName());

    // Insert the weight element.
    xmlpp::Element* pWeightElt
      = pModMolElt->add_child(eltName::weight);

    // Add the mol weight in attribute.
    double molWeight = pDefaultState->getMolWeight();
    pWeightElt->set_attribute(eltName::weight_daltonsAttr,
			      domUtils::stringify<double>(molWeight));

    // Cause all the binding sites to insert themselves.
    std::for_each(sites.begin(),
		  sites.end(),
		  insertSite(pModMolElt));

    // Run through all the modification sites, giving the name and default
    // modification of each.  This could be done with binary for_each.
    for(int modSiteNdx = 0;
	modSiteNdx < modSiteCount();
	++modSiteNdx)
      {
	// Insert element for giving the name of the modification site.
	xmlpp::Element* pModSiteElt
	  = pModMolElt->add_child(eltName::modSite);

	pModSiteElt->set_attribute(eltName::modSite_nameAttr,
				   modSiteNames[modSiteNdx]);

	// Insert element for default modification.
	xmlpp::Element* pDefaultModRefElt
	  = pModSiteElt->add_child(eltName::defaultModRef);
      
	const modification* pDefaultMod = (*pDefaultState)[modSiteNdx];
	pDefaultModRefElt->set_attribute(eltName::defaultModRef_nameAttr,
					 pDefaultMod->name);
      }

    // Run through all the registered states of this mol, displaying the
    // "allosteric" ones; i.e. all except the default state.  For each
    // (non-default) state, we only want to see the modification sites that
    // do not have their default modification.
    std::vector<siteParam> defaultSiteParams = getDefaultSiteParams();
    std::for_each(alloMap.begin(),
		  alloMap.end(),
		  insertNonDefault(pModMolElt,
				   this,
				   pDefaultState,
				   defaultSiteParams));

    return pModMolElt;
  }

  modMol*
  mustBeModMol(xmlpp::Node* pRequestingNode,
	      mol* pMol)
    throw(badModMolCastXcpt)
  {
    modMol* pModMol
      = dynamic_cast<modMol*>(pMol);

    if(! pModMol)
      throw badModMolCastXcpt(pRequestingNode,
			      pMol);

    return pModMol;
  }

  const modMol*
  mustBeModMol(xmlpp::Node* pRequestingNode,
	       const mol* pMol)
    throw(badModMolCastXcpt)
  {
    const modMol* pModMol
      = dynamic_cast<const modMol*>(pMol);

    if(! pModMol)
      throw badModMolCastXcpt(pRequestingNode,
			      pMol);

    return pModMol;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008  Walter Lawrence (Larry) Lok.
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



#include "plex/plexUnit.hh"
#include "plex/dupNodeOmniXcpt.hh"
#include "plex/noOmniForNodeXcpt.hh"

namespace plx
{
    plexUnit::
    plexUnit(mzr::moleculizer& rMoleculizer,
             mzr::mzrUnit& refMzrUnit,
             bnd::molUnit& refMolUnit,
             nmr::nmrUnit& refNmrUnit) :
        mzr::unit("plex",
                  rMoleculizer),
        rMzrUnit(refMzrUnit),
        rMolUnit(refMolUnit),
        rNmrUnit(refNmrUnit),
        recognize(*this, rNmrUnit)
    {
        // Model elements for which plex unit is responsible.
        inputCap.addModelContentName(eltName::allostericPlexes);
        inputCap.addModelContentName(eltName::allostericOmnis);

        // Explicit species for which plex unit is responsible.
        inputCap.addExplicitSpeciesContentName(eltName::plexSpecies);

        // Species streams for which plex unit is responsible.
        inputCap.addSpeciesStreamsContentName
            (eltName::plexSpeciesStream);
        inputCap.addSpeciesStreamsContentName
            (eltName::omniSpeciesStream);
    
        // Not responsible for any reaction generators.
        // Not responsible for any species.

        // The plexUnit may as well register its own omniplex Xpaths.
        const std::string slash("/");

        std::ostringstream alloOmnisXpath;
        alloOmnisXpath << mzr::eltName::model
                       << slash
                       << eltName::allostericOmnis
                       << slash
                       << eltName::allostericOmni;
        addOmniXpath(alloOmnisXpath.str());

        std::ostringstream omniSpeciesStreamXpath;
        omniSpeciesStreamXpath << mzr::eltName::streams
                               << slash
                               << mzr::eltName::speciesStreams
                               << slash
                               << eltName::omniSpeciesStream;
        addOmniXpath(omniSpeciesStreamXpath.str());
    }

    // Installs a new binding feature for the given pair of structural
    // sites. Any plex that binds the two given structural bindings
    // together will need use the bindingFeature associated to the pair
    // of structural sites to inform the family of decompositions of
    // the structural binding.
    fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> >*
    plexUnit::addBindingFeature(bnd::mzrMol* pLeftMol,
                                int leftMolSiteSpec,
                                bnd::mzrMol* pRightMol,
                                int rightMolSiteSpec)
    {
        cpx::structuralSite<bnd::mzrMol> leftSite(pLeftMol,
                                                  leftMolSiteSpec);

        cpx::structuralSite<bnd::mzrMol> rightSite(pRightMol,
                                                   rightMolSiteSpec);

        cpx::structuralBinding<bnd::mzrMol> sBinding(leftSite,
                                                     rightSite);

        fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> > 
            bindingFtr;

        // I note (13Jul05) that I don't take exception to failure of insertion.
        cpx::knownBindings<bnd::mzrMol, fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> > >::iterator
            iEntry = bindingFeatures.insert(std::make_pair(sBinding,
                                                           bindingFtr)).first;

        return &(iEntry->second);
    }

    // I expect this routine to be used in the constructor for plexFamily,
    // to construct the plexFamily's feature map for bindings.
    //
    // This routine looks for the feature "in both directions."
    fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> >*
    plexUnit::findBindingFeature(bnd::mzrMol* pLeftMol,
                                 int leftMolSiteSpec,
                                 bnd::mzrMol* pRightMol,
                                 int rightMolSiteSpec)
    {
        cpx::structuralSite<bnd::mzrMol> leftSite(pLeftMol,
                                                  leftMolSiteSpec);

        cpx::structuralSite<bnd::mzrMol> rightSite(pRightMol,
                                                   rightMolSiteSpec);

        cpx::structuralBinding<bnd::mzrMol> sBinding(leftSite,
                                                     rightSite);

        cpx::knownBindings<bnd::mzrMol, fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> > >::iterator
            iEntry = bindingFeatures.find(sBinding);

        // May have to try the sites in the opposite order.
        if(iEntry == bindingFeatures.end())
        {
            cpx::structuralBinding<bnd::mzrMol> oppBinding(rightSite,
                                                           leftSite);
            iEntry = bindingFeatures.find(oppBinding);
        }

        return iEntry == bindingFeatures.end()
            ? 0
            : &(iEntry->second);
    }

    void
    plexUnit::
    addOmniPlex(mzrOmniPlex* pOmniPlex,
                xmlpp::Node* pParentNode)
        throw(utl::xcpt)
    {
        // Many omniplexes may have the same plexFamily, so
        // we should pay no attention if this insertion fails.
        omniPlexFamilies.insert(pOmniPlex->getFamily());

        // If there is already an entry for the node, then something
        // is internally exceptional.
        if(! (nodeToOmni.insert(std::make_pair(pParentNode,
                                               pOmniPlex)).second))
            throw dupNodeOmniXcpt(pParentNode);
    }

    mzrOmniPlex*
    plexUnit::
    getOmniForNode(xmlpp::Node* pParentNode) const
    {
        std::map<xmlpp::Node*, mzrOmniPlex*>::const_iterator
            iEntry
            = nodeToOmni.find(pParentNode);

        return (nodeToOmni.end() == iEntry)
            ? 0
            : iEntry->second;
    }

    mzrOmniPlex*
    plexUnit::
    mustGetOmniForNode(xmlpp::Node* pParentNode) const
        throw(utl::xcpt)
    {
        mzrOmniPlex* pOmni
            = getOmniForNode(pParentNode);

        if(! pOmni) throw noOmniForNodeXcpt(pParentNode);

        return pOmni;
    }

    // For dumping state in XML.
    void
    plexUnit::
    insertStateElts(xmlpp::Element* pRootElt)
        throw(std::exception)
    {
        // Get the model element.
        xmlpp::Element* pModelElt
            = utl::dom::mustGetUniqueChild(pRootElt,
                                           mzr::eltName::model);

        // Ensure that the tagged-species element is here.
        xmlpp::Element* pTaggedSpeciesElt
            = utl::dom::mustGetUniqueChild(pModelElt,
                                           mzr::eltName::taggedSpecies);

        // Insert tagged-plex-species nodes.
        recognize.insertSpecies(pTaggedSpeciesElt,
                                rMzrUnit.getMolarFactor().getFactor());
    }

    unsigned int
    plexUnit::
    getNumberPlexSpecies() const
    {
        throw "Help from plexUnit::getNumberPlexSpecies";
        return 0;

    }
    

    mzrPlexSpecies*
    plexUnit::constructNewPlexSpeciesFromComplexOutputState(nmr::ComplexOutputStateCref aCOS)
    {
        plx::mzrPlex* pMzrPlex = new plx::mzrPlex;

        for( std::vector<std::string>::const_iterator iter = aCOS.theMolTokens.begin();
             iter != aCOS.theMolTokens.end();
             ++iter)
        {
            try
            {
                bnd::mzrMol* ptrMzrMol = rMolUnit.mustFindMol( *iter );
                pMzrPlex->mols.push_back( ptrMzrMol );

                std::cout << "########  BEGIN ##############################" << std::endl;
                std::cout << "Mol: " << ptrMzrMol->getName() << std::endl;
                for( std::vector<bnd::mzrBndSite>::const_iterator iter = ((std::vector<bnd::mzrBndSite>*) ptrMzrMol)->begin();
                     iter != ((std::vector<bnd::mzrBndSite>*) ptrMzrMol)->end();
                     ++iter)
                {
                    std::cout << "\t" << iter->getName() << std::endl;
                }

                std::cout << "#######   END  ##############################" << std::endl;
            }

            catch( nmr::MissingNameEncoderXcpt e)
            {
                e.wailAndBail();
            }
            catch(utl::xcpt e)
            {
                e.wailAndBail();
            }
        }

        for(std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > >::const_iterator iter = aCOS.theBindingTokens.begin();
            iter != aCOS.theBindingTokens.end();
            ++iter)
        {
            int first_molIndex, first_bndIndex, second_molIndex, second_bndIndex;

            utl::from_string( first_molIndex, iter->first.first);
            utl::from_string( first_bndIndex, iter->first.second);
            utl::from_string( second_molIndex, iter->second.first);
            utl::from_string( second_bndIndex, iter->second.second);

            cpx::siteSpec firstBinding( first_molIndex, first_bndIndex);
            cpx::siteSpec secondBinding( second_molIndex, second_bndIndex);

            pMzrPlex->bindings.push_back( cpx::binding( firstBinding, secondBinding) );

//             std::cout << "Creating binding from " << pMzrPlex->mols[first_molIndex]->getName() << ", " 
//                       << (*pMzrPlex->mols[first_molIndex])[first_bndIndex].getName() 
//                       << " to " 
//                       << pMzrPlex->mols[second_molIndex]->getName() << ", " 
//                       << (*pMzrPlex->mols[second_molIndex])[second_bndIndex].getName()  << std::endl;
        }

        cpx::plexIso resultToResultParadigm;

        // Recognize to create a mzrPlexFamily.
        // This takes ownership of the pMzrPlex.
	plx::mzrPlexFamily* pProductFamily
            = recognize(*pMzrPlex, resultToResultParadigm);

        std::vector<cpx::molParam> theDefaultParams = pProductFamily->makeDefaultMolParams();

        // TODO -- lots of todos
        for (std::vector<std::pair< std::string, std::pair<std::string, std::string> > >::const_iterator iter = aCOS.theModificationTokens.begin();
             iter != aCOS.theModificationTokens.end();
             ++iter)
        {
            int molIndex, modificationIndex;
            std::string modificationState( iter->second.second);
            utl::from_string(molIndex, iter->first);
            utl::from_string(modificationIndex, iter->second.first);

            // Make sure this is correct, as opposed to the backwards map.
            int mappedMolIndex = resultToResultParadigm.forward.applyToMolSpec( molIndex );

            // the mod with index 'mappedMolIndex' should have modification 'modificationState'
            // applied to itself at modification site numbered 'modificationIndex'.


            
        }

        return NULL;

    }
    
}
        

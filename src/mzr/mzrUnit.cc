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

#include "utl/stdIncludes.hh"

#include "mzr/mzrUnit.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/moleculizer.hh"
#include "mzr/dupSpeciesNameXcpt.hh"
#include "mzr/unkSpeciesXcpt.hh"
#include "mzr/dupDumpableNameXcpt.hh"
#include "mzr/unkDumpableXcpt.hh"
#include "plex/plexUnit.hh"

namespace mzr
{


    void
    mzrUnit::
    mustAddSpecies(const std::string& rSpeciesName,
                   mzrSpecies* pSpecies,
                   xmlpp::Node* pRequestingNode)
        throw(utl::xcpt)
    {
        if(! addSpecies(rSpeciesName,
                        pSpecies))
            throw dupSpeciesNameXcpt(rSpeciesName,
                                     pRequestingNode);
    }

    mzrSpecies*
    mzrUnit::
    mustFindSpecies(const std::string& rSpeciesName,
                    xmlpp::Node* pRequestingNode) const
        throw(utl::xcpt)
    {
        mzrSpecies* pSpecies = findSpecies(rSpeciesName);

        if(! pSpecies)
            throw unkSpeciesXcpt(rSpeciesName,
                                 pRequestingNode);
        return pSpecies;
    }

    void
    mzrUnit::
    mustAddDumpable(fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable,
                    xmlpp::Node* pRequestingNode)
        throw(utl::xcpt)
    {
        if(! addDumpable(pDumpable))
            throw dupDumpableNameXcpt(pDumpable->getName(),
                                      pRequestingNode);
    }
  
    fnd::dumpable<fnd::basicDumpable::dumpArg>*
    mzrUnit::
    mustFindDumpable(const std::string& rDumpableName,
                     xmlpp::Node* pRequestingNode) const
        throw(utl::xcpt)
    {
        fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable
            = findDumpable(rDumpableName);

        if(! pDumpable)
            throw unkDumpableXcpt(rDumpableName,
                                  pRequestingNode);

        return pDumpable;
    }



    mzrUnit::
    mzrUnit(moleculizer& rMoleculizer) :
        unit("mzr",
             rMoleculizer),
        generateDepth(0)
    {
        // Initialize the vector of global state variables, with which each
        // reaction is initialized.  At this point, there is only one
        // global state variable, the volume.
        globalVars.push_back(&theMolarFactor);
    
        // Model elements whose contents are parsed by some unit
        // or another, as determined by moleculizer::parseDomInput.
        inputCap.addModelContentName(eltName::reactionGens);
        inputCap.addModelContentName(eltName::explicitSpecies);

        // Model elements that this unit actually processes.
        inputCap.addModelContentName(eltName::explicitReactions);
        inputCap.addModelContentName(eltName::volume);

        // This unit is not responsible for any reaction generators
        // or species streams.
    }


    void
    mzrUnit::prepareToRun(xmlpp::Element* pRootElt,
                          xmlpp::Element* pModelElt) throw(std::exception)
    {}

    void
    mzrUnit::prepareToContinue(xmlpp::Element* pRootElt,
                               xmlpp::Element* pModelElt,
                               xmlpp::Element* pStreamsElt,
                               std::map<std::string, std::string>& rTagToName,
                               xmlpp::Element* pTaggedSpeciesElement)
        throw(std::exception)
    {}


}

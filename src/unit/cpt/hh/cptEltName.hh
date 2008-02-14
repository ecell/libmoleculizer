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

#ifndef CPTELTNAME_H
#define CPTELTNAME_H

#include <string>

namespace cpt
{
  namespace eltName
  {
    // Names drawn from mzr::mzrEltName.  These should be substantially
    // the same for cpt when cpt is filled out.

    extern const std::string cptInput;
    extern const std::string model;
    extern const std::string streams;
    extern const std::string events;

    extern const std::string explicitSpecies;

    extern const std::string speciesStreams;

    extern const std::string dumpStreams;
    extern const std::string dumpStream;
    extern const std::string dumpStream_dumpPeriodAttr;
    extern const std::string targetFile;
    extern const std::string targetFile_fileNameAttr;

    extern const std::string speciesStreamRef;
    extern const std::string speciesStreamRef_nameAttr;
    extern const std::string speciesStreamRef_totalPopsAttr;
    extern const std::string speciesStreamRef_yesTotalPops;
    extern const std::string speciesStreamRef_noTotalPops;

    extern const std::string speciesRef;
    extern const std::string speciesRef_nameAttr;
    extern const std::string speciesRef_totalPopsAttr;
    extern const std::string speciesRef_yesTotalPops;
    extern const std::string speciesRef_noTotalPops;

    extern const std::string statStreamRef;
    extern const std::string statStreamRef_nameAttr;
    extern const std::string statStream_simTime;
    extern const std::string statStream_clockTime;
    extern const std::string statStream_clockSeconds;
    extern const std::string statStream_activationCount;
    extern const std::string statStream_speciesCount;
    extern const std::string statStream_reactionCount;
    extern const std::string statStream_reactionEventCount;
    extern const std::string statStream_volume;

    extern const std::string reactionGens;

    extern const std::string explicitReactions;
    extern const std::string reaction;
    extern const std::string substrateSpeciesRef;
    extern const std::string substrateSpeciesRef_nameAttr;
    extern const std::string substrateSpeciesRef_multAttr;
    extern const std::string productSpeciesRef;
    extern const std::string productSpeciesRef_nameAttr;
    extern const std::string productSpeciesRef_multAttr;
    extern const std::string rate;
    extern const std::string rate_valueAttr;

    extern const std::string volume;
    extern const std::string volume_litersAttr;

    extern const std::string createEvent;
    extern const std::string createEvent_timeAttr;
    extern const std::string population;
    extern const std::string population_countAttr;
    extern const std::string dumpStreamRef;
    extern const std::string dumpStreamRef_nameAttr;
    extern const std::string dumpStateEvent;
    extern const std::string dumpStateEvent_timeAttr;

    extern const std::string noReactEvent;
    extern const std::string noReactEvent_timeAttr;

    extern const std::string stopEvent;
    extern const std::string stopEvent_timeAttr;

    // Names from the first, experimental input format for
    // the compartmental simulator.

    extern const std::string compartmentModel;

    extern const std::string compartmentGraph;
    extern const std::string compartments;
    extern const std::string compartment;
    extern const std::string compartment_nameAttr;
    extern const std::string compartment_volumeAttr;
    extern const std::string boundaries;
    extern const std::string boundary;
    extern const std::string boundary_leftCompartmentAttr;
    extern const std::string boundary_rightCompartmentAttr;
    extern const std::string boundary_areaAttr;
    extern const std::string boundary_thicknessAttr;

    extern const std::string dumpCompartments;
    extern const std::string compartmentRef;
    extern const std::string compartmentRef_nameAttr;

    extern const std::string explicitGlobalSpecies;
    extern const std::string globalSpecies;
    extern const std::string globalSpecies_nameAttr;
    extern const std::string defaultDiffusionRate;
    extern const std::string defaultDiffusionRate_valueAttr;
    extern const std::string boundaryRate;
    extern const std::string boundaryRate_valueAttr;
    
    extern const std::string defaultCompartmentPop;
    extern const std::string defaultCompartmentPop_countAttr;
    extern const std::string compartmentPop;
    extern const std::string compartmentPop_compartmentNameAttr;
    extern const std::string compartmentPop_countAttr;

    extern const std::string explicitGlobalReactions;
    extern const std::string globalReaction;
    extern const std::string reactantGlobalSpeciesRef;
    extern const std::string reactantGlobalSpeciesRef_nameAttr;
    extern const std::string reactantGlobalSpeciesRef_multAttr;
    extern const std::string productGlobalSpeciesRef;
    extern const std::string productGlobalSpeciesRef_nameAttr;
    extern const std::string productGlobalSpeciesRef_multAttr;

    extern const std::string runParams;
    extern const std::string runParams_cycleTimeAttr;
    extern const std::string runParams_reactionsPerCycleAttr;
    extern const std::string runParams_epsilonAttr;
    extern const std::string runParams_runTimeAttr;

    // Elements pertaining to output.

    // Names for dumping scientific notation for double-valued parameters,
    // intentded for use in SBML, etc.
    extern const std::string sciNote;
    extern const std::string sciNote_fractionAttr;
    extern const std::string sciNote_exponentAttr;

    extern const std::string moleculizerState;
    extern const std::string unitsStates;

    extern const std::string explicitSpeciesTags;
    extern const std::string explicitSpeciesTag;
    extern const std::string explicitSpeciesTag_nameAttr;
    extern const std::string explicitSpeciesTag_tagAttr;

    extern const std::string taggedSpecies;

    extern const std::string tagReactions;
    extern const std::string tagReaction;
    extern const std::string taggedSubstrate;
    extern const std::string taggedSubstrate_tagAttr;
    extern const std::string taggedSubstrate_multAttr;
    extern const std::string taggedProduct;
    extern const std::string taggedProduct_tagAttr;
    extern const std::string taggedProduct_multAttr;

    extern const std::string time;
    extern const std::string time_secondsAttr;

    extern const std::string taggedSpeciesStreams;
    extern const std::string taggedSpeciesStream;
    extern const std::string taggedSpeciesStream_nameAttr;
    extern const std::string taggedSpeciesRef;
    extern const std::string taggedSpeciesRef_tagAttr;

    extern const std::string taggedDumpStreams;
    extern const std::string taggedDumpStream;
    extern const std::string taggedDumpStream_fileNameAttr;
    extern const std::string taggedDumpStream_dumpPeriodAttr;
    extern const std::string taggedSpeciesStreamRef;
    extern const std::string taggedSpeciesStreamRef_nameAttr;

  }
}

#endif // CPTELTNAME_H

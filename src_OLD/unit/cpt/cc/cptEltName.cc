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

#include "cpt/cptEltName.hh"

namespace cpt
{
  namespace eltName
  {
    // Names drawn from mzr::mzrEltName.  These should be substantially
    // the same for cpt when cpt is filled out.

    const std::string cptInput("cpt-input");
    const std::string model("model");
    const std::string streams("streams");
    const std::string events("events");

    const std::string explicitSpecies("explicit-species");

    const std::string speciesStreams("species-streams");

    const std::string dumpStreams("dump-streams");
    const std::string dumpStream("dump-stream");
    const std::string dumpStream_dumpPeriodAttr("dump-period");
    const std::string targetFile("target-file");
    const std::string targetFile_fileNameAttr("file-name");

    const std::string speciesStreamRef("species-stream-ref");
    const std::string speciesStreamRef_nameAttr("name");
    const std::string speciesStreamRef_totalPopsAttr("total-pops");
    const std::string speciesStreamRef_yesTotalPops("yes");
    const std::string speciesStreamRef_noTotalPops("no");
      
    const std::string speciesRef("species-ref");
    const std::string speciesRef_nameAttr("name");
    const std::string speciesRef_totalPopsAttr("total-pops");
    const std::string speciesRef_yesTotalPops("yes");
    const std::string speciesRef_noTotalPops("no");

    const std::string statStreamRef("stat-stream-ref");
    const std::string statStreamRef_nameAttr("name");
    const std::string statStream_simTime("sim-time");
    const std::string statStream_clockTime("clock-time");
    const std::string statStream_clockSeconds("clock-seconds");
    const std::string statStream_activationCount("activation-count");
    const std::string statStream_speciesCount("species-count");
    const std::string statStream_reactionCount("reaction-count");
    const std::string statStream_reactionEventCount("reaction-event-count");
    const std::string statStream_volume("volume");
//     const std::string statStream_firstReactionQueueSize("first-reaction-queue-size");
//     const std::string statStream_nextReactionQueueSize("next-reaction-queue-size");
    

    const std::string reactionGens("reaction-gens");

    const std::string explicitReactions("explicit-reactions");
    const std::string reaction("reaction");
    const std::string substrateSpeciesRef("substrate-species-ref");
    const std::string substrateSpeciesRef_nameAttr("name");
    const std::string substrateSpeciesRef_multAttr("multiplicity");
    const std::string productSpeciesRef("product-species-ref");
    const std::string productSpeciesRef_nameAttr("name");
    const std::string productSpeciesRef_multAttr("multiplicity");

    const std::string volume("volume");
    const std::string volume_litersAttr("liters");

    const std::string createEvent("create-event");
    const std::string createEvent_timeAttr("time");
    const std::string population("population");
    const std::string population_countAttr("count");
    const std::string dumpStreamRef("dump-stream-ref");
    const std::string dumpStreamRef_nameAttr("name");
    const std::string dumpStateEvent("dump-state-event");
    const std::string dumpStateEvent_timeAttr("time");

    const std::string noReactEvent("no-react-event");
    const std::string noReactEvent_timeAttr("time");

    const std::string stopEvent("stop-event");
    const std::string stopEvent_timeAttr("time");

    // Names from the first, experimental input format for
    // the compartmental simulator.

    const std::string compartmentModel("compartment-model");

    const std::string compartmentGraph("compartment-graph");
    const std::string compartments("compartments");
    const std::string compartment("compartment");
    const std::string compartment_nameAttr("name");
    const std::string compartment_volumeAttr("volume");
    const std::string boundaries("boundaries");
    const std::string boundary("boundary");
    const std::string boundary_leftCompartmentAttr("left-compartment");
    const std::string boundary_rightCompartmentAttr("right-compartment");
    const std::string boundary_areaAttr("area");
    const std::string boundary_thicknessAttr("thickness");

    const std::string dumpCompartments("dump-compartments");
    const std::string compartmentRef("compartment-ref");
    const std::string compartmentRef_nameAttr("name");

    const std::string explicitGlobalSpecies("explicit-global-species");
    const std::string globalSpecies("global-species");
    const std::string globalSpecies_nameAttr("name");
    const std::string defaultDiffusionRate("default-diffusion-rate");
    const std::string defaultDiffusionRate_valueAttr("value");
    const std::string boundaryRate("boundary-rate");
    const std::string boundaryRate_valueAttr("value");
    const std::string defaultCompartmentPop("default-compartment-pop");
    const std::string defaultCompartmentPop_countAttr("count");
    const std::string compartmentPop("compartment-pop");
    const std::string compartmentPop_compartmentNameAttr("compartment-name");
    const std::string compartmentPop_countAttr("count");

    const std::string explicitGlobalReactions("explicit-global-reactions");
    const std::string globalReaction("global-reaction");
    const std::string reactantGlobalSpeciesRef("reactant-global-species-ref");
    const std::string reactantGlobalSpeciesRef_nameAttr("name");
    const std::string reactantGlobalSpeciesRef_multAttr("multiplicity");
    const std::string productGlobalSpeciesRef("product-global-species-ref");
    const std::string productGlobalSpeciesRef_nameAttr("name");
    const std::string productGlobalSpeciesRef_multAttr("multiplicity");
    const std::string rate("rate");
    const std::string rate_valueAttr("value");

    const std::string runParams("run-params");
    const std::string runParams_cycleTimeAttr("cycle-time");
    const std::string runParams_reactionsPerCycleAttr("reactions-per-cycle");
    const std::string runParams_epsilonAttr("epsilon");
    const std::string runParams_runTimeAttr("run-time");


    // Elements pertaining to output.

    // Names for dumping scientific notation for double-valued parameters,
    // intentded for use in SBML, etc.
    const std::string sciNote("sci-note");
    const std::string sciNote_fractionAttr("fraction");
    const std::string sciNote_exponentAttr("exponent");

    const std::string moleculizerState("moleculizer-state");
    const std::string unitsStates("units-states");

    const std::string explicitSpeciesTags("explicit-species-tags");
    const std::string explicitSpeciesTag("explicit-species-tag");
    const std::string explicitSpeciesTag_nameAttr("name");
    const std::string explicitSpeciesTag_tagAttr("tag");

    const std::string taggedSpecies("tagged-species");

    const std::string tagReactions("tag-reactions");
    const std::string tagReaction("tag-reaction");
    const std::string taggedSubstrate("tagged-substrate");
    const std::string taggedSubstrate_tagAttr("tag");
    const std::string taggedSubstrate_multAttr("multiplicity");
    const std::string taggedProduct("tagged-product");
    const std::string taggedProduct_tagAttr("tag");
    const std::string taggedProduct_multAttr("multiplicity");

    const std::string time("time");
    const std::string time_secondsAttr("seconds");

    const std::string taggedSpeciesStreams("tagged-species-streams");
    const std::string taggedSpeciesStream("tagged-species-stream");
    const std::string taggedSpeciesStream_nameAttr("name");
    const std::string taggedSpeciesRef("tagged-species-ref");
    const std::string taggedSpeciesRef_tagAttr("tag");

    const std::string taggedDumpStreams("tagged-dump-streams");
    const std::string taggedDumpStream("tagged-dump-stream");
    const std::string taggedDumpStream_fileNameAttr("file-name");
    const std::string taggedDumpStream_dumpPeriodAttr("dump-period");
    const std::string taggedSpeciesStreamRef("tagged-species-stream-ref");
    const std::string taggedSpeciesStreamRef_nameAttr("name");
  }
}

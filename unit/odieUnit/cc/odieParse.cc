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

#include "odie/odieApp.hh"
#include "odie/odieEltName.hh"
#include "odie/odieParse.hh"

namespace odie
{
  class parseSpeciesEntry :
    public std::unary_function<xmlpp::Node*, int>
  {
    odieSpeciesCatalog& rCatalog;
    double* pConc;
    
  public:
    parseSpeciesEntry(odieSpeciesCatalog& rOdieSpeciesCatalog,
		      double* pConcentrations) :
      rCatalog(rOdieSpeciesCatalog),
      pConc(pConcentrations)
    {}
  
    void
    operator()(xmlpp::Node* pSpeciesEntryNode) throw(std::exception)
    {
      xmlpp::Element* pSpeciesEntryElt
	= domUtils::mustBeElementPtr(pSpeciesEntryNode);

      int speciesNdx = rCatalog.size();

      // Get the name and population of the species.
      std::string speciesName
	= domUtils::mustGetAttrString(pSpeciesEntryElt,
				      odieEltName::speciesEntry_nameAttr);

      // Add the species's concentration to the buffer.
      pConc[speciesNdx]
	= domUtils::mustGetAttrDouble(pSpeciesEntryElt,
				      odieEltName::speciesEntry_concAttr);

      // Add an entry to the catalog so we can find the species's
      // index from its name later in the parse.
      rCatalog.mustAddUnique(pSpeciesEntryNode,
			     speciesName,
			     speciesNdx);
    }
  };

  class parseReactionSubstrateNode :
    public std::unary_function<xmlpp::Node*, std::map<int, int>::value_type>
  {
    odieSpeciesCatalog& rSpeciesCatalog;
    
  public:
    parseReactionSubstrateNode(odieSpeciesCatalog& rOdieSpeciesCatalog) :
      rSpeciesCatalog(rOdieSpeciesCatalog)
    {}

    result_type
    operator()(xmlpp::Node* pReactionSubstrateNode) throw(std::exception)
    {
      xmlpp::Element* pReactionSubstrateElt
	= domUtils::mustBeElementPtr(pReactionSubstrateNode);

      std::string substrateName
	= domUtils::mustGetAttrString(pReactionSubstrateElt,
				      odieEltName::reactionSubstrate_nameAttr);
      int speciesNdx
	= rSpeciesCatalog.mustFindSpecies(pReactionSubstrateElt,
					  substrateName);
      int multiplicity
	= domUtils::mustGetAttrInt(pReactionSubstrateElt,
				   odieEltName::reactionSubstrate_multAttr);

      return std::make_pair(speciesNdx,
			    multiplicity);
    }
  };

  class parseReactionProductNode :
    public std::unary_function<xmlpp::Node*, std::map<int, int>::value_type>
  {
    odieSpeciesCatalog& rSpeciesCatalog;

  public:
    parseReactionProductNode(odieSpeciesCatalog& rOdieSpeciesCatalog) :
      rSpeciesCatalog(rOdieSpeciesCatalog)
    {}

    result_type
    operator()(xmlpp::Node* pReactionProductNode) throw(std::exception)
    {
      xmlpp::Element* pReactionProductElt
	= domUtils::mustBeElementPtr(pReactionProductNode);

      std::string productSpeciesName
	= domUtils::mustGetAttrString(pReactionProductElt,
				      odieEltName::reactionProduct_nameAttr);
      int speciesNdx
	= rSpeciesCatalog.mustFindSpecies(pReactionProductElt,
					  productSpeciesName);
      int multiplicity
	= domUtils::mustGetAttrInt(pReactionProductElt,
				   odieEltName::reactionProduct_multAttr);

      return std::make_pair(speciesNdx,
			    multiplicity);
    }
  };

  class parseReaction :
    public std::unary_function<xmlpp::Node*, parserReaction>
  {
    odieSpeciesCatalog& rSpeciesCatalog;
  public:
    parseReaction(odieSpeciesCatalog& rOdieSpeciesCatalog) :
      rSpeciesCatalog(rOdieSpeciesCatalog)
    {}
  
    parserReaction
    operator()(xmlpp::Node* pReactionNode) throw(std::exception)
    {
      parserReaction reaction;

      xmlpp::Element* pReactionElt
	= domUtils::mustBeElementPtr(pReactionNode);

      // Parse the substrates and their multiplicities.
      xmlpp::Node::NodeList substrateNodes
	= pReactionElt->get_children(odieEltName::reactionSubstrate);
      std::transform(substrateNodes.begin(),
		     substrateNodes.end(),
		     std::inserter(reaction.substrates,
				   reaction.substrates.begin()),
		     parseReactionSubstrateNode(rSpeciesCatalog));

      // Parse the products (species and multiplicities just like substrates).
      xmlpp::Node::NodeList productNodes
	= pReactionElt->get_children(odieEltName::reactionProduct);
      std::transform(productNodes.begin(),
		     productNodes.end(),
		     std::inserter(reaction.products,
				   reaction.products.begin()),
		     parseReactionProductNode(rSpeciesCatalog));

      // Parse the reaction rate.
      xmlpp::Element* pRateElt
	= domUtils::mustGetUniqueChild(pReactionElt,
				       odieEltName::rate);
      reaction.rate
	= domUtils::mustGetAttrDouble(pRateElt,
				      odieEltName::rate_valueAttr);

      return reaction;
    }
  };

  class lookUpDumpableSpecies :
    public std::unary_function<xmlpp::Node*, int>
  {
    odieSpeciesCatalog& rSpeciesCatalog;
  public:
    lookUpDumpableSpecies(odieSpeciesCatalog& rOdieSpeciesCatalog) :
      rSpeciesCatalog(rOdieSpeciesCatalog)
    {}

    int
    operator()(xmlpp::Node* pSpeciesRefNode) throw(std::exception)
    {
      xmlpp::Element* pSpeciesRefElt
	= domUtils::mustBeElementPtr(pSpeciesRefNode);

      std::string speciesName
	= domUtils::mustGetAttrString(pSpeciesRefElt,
				      odieEltName::speciesRef_nameAttr);

      return rSpeciesCatalog.mustFindSpecies(pSpeciesRefElt,
					     speciesName);
    }
  };

  // Parses a "dumpable" node, which contains the name of the dumpable,
  // destined to be a column header in a dump file, along with a list
  // of one or more species.
  class parseDumpableNode :
    public std::unary_function<xmlpp::Node*, odieDumpableCatalog::value_type>
  {
    odieSpeciesCatalog& rSpeciesCatalog;
  public:
    parseDumpableNode(odieSpeciesCatalog& rOdieSpeciesCatalog) :
      rSpeciesCatalog(rOdieSpeciesCatalog)
    {}

    result_type
    operator()(xmlpp::Node* pDumpableNode) throw(std::exception)
    {
      odieDumpable dumpable;
    
      xmlpp::Element* pDumpableElt
	= domUtils::mustBeElementPtr(pDumpableNode);

      dumpable.name
	= domUtils::mustGetAttrString(pDumpableElt,
				      odieEltName::dumpable_nameAttr);

      xmlpp::Node::NodeList speciesRefNodes
	= pDumpableElt->get_children(odieEltName::speciesRef);
      std::transform(speciesRefNodes.begin(),
		     speciesRefNodes.end(),
		     back_inserter(dumpable),
		     lookUpDumpableSpecies(rSpeciesCatalog));

      return std::make_pair(dumpable.name,
			    dumpable);
    }
  };

  class parseDumpableRefNode :
    public std::unary_function<xmlpp::Node*, odieDumpable>
  {
    const odieDumpableCatalog& rDumpableCatalog;
  public:
    parseDumpableRefNode(const odieDumpableCatalog& rOdieDumpableCatalog) :
      rDumpableCatalog(rOdieDumpableCatalog)
    {}
  
    const odieDumpable&
    operator()(xmlpp::Node* pDumpableRefNode) throw(std::exception)
    {
      xmlpp::Element* pDumpableRefElt
	= domUtils::mustBeElementPtr(pDumpableRefNode);

      std::string dumpableName
	= domUtils::mustGetAttrString(pDumpableRefElt,
				      odieEltName::dumpableRef_nameAttr);

      return rDumpableCatalog.mustFindDumpable(pDumpableRefElt,
					       dumpableName);
    }
  };

  class parseDumpStreamNode :
    public std::unary_function<xmlpp::Node*, dumpStream>
  {
    const odieDumpableCatalog& rDumpableCatalog;
  public:
    parseDumpStreamNode(const odieDumpableCatalog& rOdieDumpableCatalog) :
      rDumpableCatalog(rOdieDumpableCatalog)
    {}

    dumpStream
    operator()(xmlpp::Node* pDumpStreamNode) throw(std::exception)
    {
      dumpStream result;
    
      xmlpp::Element* pDumpStreamElt
	= domUtils::mustBeElementPtr(pDumpStreamNode);

      // Parse the file name, which seems to be given by the "name"
      // attribute of the dumpStream.
      std::string fileName
	= domUtils::mustGetAttrString(pDumpStreamElt,
				      odieEltName::dumpStream_nameAttr);

      // Parse the dumpables that are supposed to appear as columns
      // in the output file.
      xmlpp::Node::NodeList dumpableRefNodes
	= pDumpStreamElt->get_children(odieEltName::dumpableRef);

      std::transform(dumpableRefNodes.begin(),
		     dumpableRefNodes.end(),
		     std::back_inserter(result),
		     parseDumpableRefNode(rDumpableCatalog));

      result.mustSetStream(pDumpStreamElt,
			   fileName);
      return result;
    }
  };

  odieApp::odieApp(int argc,
		   char** argv,
		   xmlpp::Document* pDoc) throw(std::exception) :
    pConcentrations(0),
    now(0.0),
    stopTime(0.0),
    epsilonAbs(1.0),
    epsilonRel(1.0),
    stateScalar(1.0),
    derivScalar(1.0)
  {
    xmlpp::Element* pRootElt
      = pDoc->get_root_node();

    xmlpp::Element* pSpeciesElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::species);

    // Parse the species names and concentrations.

    // Get the nodes that yield species.
    xmlpp::Node::NodeList speciesEntryNodes
      = pSpeciesElt->get_children(odieEltName::speciesEntry);

    // Allocate room for the species's concentrations.
    // We can't push them back directly; an alternative would
    // be to accumualate them into a vector, then allocate
    // and fill the pConcentrations buffer.
    pConcentrations = new double[speciesEntryNodes.size()];

    // Resize the derivative.
    derivative = rk4util::polymap<double>(speciesEntryNodes.size());

    odieSpeciesCatalog speciesCatalog;

    std::for_each(speciesEntryNodes.begin(),
		  speciesEntryNodes.end(),
		  parseSpeciesEntry(speciesCatalog,
				    pConcentrations));

    // Parse the reactions, using the species catalog
    // to look up substrate and product species.
    xmlpp::Element* pReactionsElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::reactions);
    
    xmlpp::Node::NodeList reactionNodes
      = pReactionsElt->get_children(odieEltName::reaction);

    std::vector<parserReaction> reactions;
    std::transform(reactionNodes.begin(),
		   reactionNodes.end(),
		   std::back_inserter(reactions),
		   parseReaction(speciesCatalog));

    // Parse the dumpables.  These are used to cook up dumpStreams, which
    // are all that the program knows.
    xmlpp::Element* pDumpablesElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::dumpables);

    xmlpp::Node::NodeList dumpableNodes
      = pDumpablesElt->get_children(odieEltName::dumpable);
    odieDumpableCatalog dumpables;
    std::transform(dumpableNodes.begin(),
		   dumpableNodes.end(),
		   std::inserter(dumpables,
				 dumpables.begin()),
		   parseDumpableNode(speciesCatalog));

    // Parse the dump-streams.
    xmlpp::Element* pDumpStreamsElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::dumpStreams);

    xmlpp::Node::NodeList dumpStreamNodes
      = pDumpStreamsElt->get_children(odieEltName::dumpStream);

    // The name of a dumpStream is really a file name.  In the old parser code,
    // the file was opened when the name of the stream was read.
    std::transform(dumpStreamNodes.begin(),
		   dumpStreamNodes.end(),
		   std::back_inserter(dumpStreams),
		   parseDumpStreamNode(dumpables));

    // Parse time to stop simulation.
    xmlpp::Element* pStopTimeElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::stopTime);
    stopTime
      = domUtils::mustGetAttrDouble(pStopTimeElt,
				    odieEltName::stopTime_secondsAttr);
    

    // Parse control parameters for adaptive step size.
    xmlpp::Element* pEpsilonAbsElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::epsilonAbs);
    epsilonAbs
      = domUtils::mustGetAttrDouble(pEpsilonAbsElt,
				    odieEltName::epsilonAbs_valueAttr);

    xmlpp::Element* pEpsilonRelElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::epsilonRel);

    epsilonRel
      = domUtils::mustGetAttrDouble(pEpsilonRelElt,
				    odieEltName::epsilonRel_valueAttr);

    xmlpp::Element* pStateScalarElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::stateScalar);

    stateScalar
      = domUtils::mustGetAttrDouble(pStateScalarElt,
				    odieEltName::stateScalar_valueAttr);

    xmlpp::Element* pDerivScalarElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     odieEltName::derivScalar);

    derivScalar
      = domUtils::mustGetAttrDouble(pDerivScalarElt,
				    odieEltName::derivScalar_valueAttr);

    // Create the differential equation from the reactions.
    makeDerivatives(reactions);
  }
}

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

#include "utl/dom.hh"
#include "fnd/physConst.hh"
#include "rk4/tauEltName.hh"
#include "rk4/dump.hh"
#include "rk4/tauApp.hh"
#include "rk4/tauParse.hh"
#include "rk4/unkDumpableXcpt.hh"
#include "rk4/unkSpeciesXcpt.hh"

namespace rk4tau
{
  class duplicateSpeciesEntryXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingNode,
	  const std::string duplicateSpeciesName)
    {
      std::ostringstream msgStream;
      msgStream << utl::dom::xcpt::mkMsg(pOffendingNode)
		<< "Duplicate entry for species `"
		<< duplicateSpeciesName
		<< "'.";
      return msgStream.str();
    }
  public:
    duplicateSpeciesEntryXcpt(xmlpp::Node* pOffendingNode,
			      const std::string duplicateSpeciesName) :
      utl::xcpt(mkMsg(pOffendingNode,
		      duplicateSpeciesName))
      {}
  };
	
  class tauDumpableCatalog : public std::map<std::string,  tauDumpable>
  {
  public:
    // Consider returning a reference here?
    const tauDumpable&
    mustFindDumpable(xmlpp::Node* pRequestingNode,
		     const std::string& rDumpableName) const
    {
      const_iterator iEntry = find(rDumpableName);
      if(iEntry == end()) throw unkDumpableXcpt(rDumpableName,
						pRequestingNode);
      return iEntry->second;
    }
  };

  class tauSpeciesCatalog : public std::map<std::string, int>
  {
  public:
    void
    mustAddUnique(xmlpp::Node* pRequestingNode,
		  const std::string& rSpeciesName,
		  int speciesNdx)
      throw(duplicateSpeciesEntryXcpt)
    {
      if(! insert(std::make_pair(rSpeciesName,
				 speciesNdx)).second)
	throw duplicateSpeciesEntryXcpt(pRequestingNode,
					rSpeciesName);
    }

    int
    mustFindSpecies(xmlpp::Node* pRequestingNode,
		    const std::string& rSpeciesName) const
      throw(unkSpeciesXcpt)
    {
      const_iterator iEntry = find(rSpeciesName);
      if(iEntry == end()) throw unkSpeciesXcpt(rSpeciesName,
					       pRequestingNode);
      return iEntry->second;
    }
  };

  class parseSpeciesEntry :
    public std::unary_function<xmlpp::Node*, int>
  {
    tauSpeciesCatalog& rSpeciesCatalog;
    std::vector<int>& rPops;
    
  public:
    parseSpeciesEntry(tauSpeciesCatalog& rTauSpeciesCatalog,
		      std::vector<int>& rPopulations) :
      rSpeciesCatalog(rTauSpeciesCatalog),
      rPops(rPopulations)
    {}
  
    void
    operator()(xmlpp::Node* pSpeciesEntryNode) throw(std::exception)
    {
      xmlpp::Element* pSpeciesEntryElt
	= utl::dom::mustBeElementPtr(pSpeciesEntryNode);

      // Get the name and population of the species.
      std::string speciesName
	= utl::dom::mustGetAttrString(pSpeciesEntryElt,
				      tauEltName::speciesEntry_nameAttr);
      int population
	= utl::dom::mustGetAttrInt(pSpeciesEntryElt,
				   tauEltName::speciesEntry_popAttr);

      int speciesNdx = rPops.size();

      // Add an entry to the catalog so we can find the species's
      // index from its name later in the parse.
      rSpeciesCatalog.mustAddUnique(pSpeciesEntryNode,
				    speciesName,
				    speciesNdx);

      rPops.push_back(population);
    }
  };

  class parseReactionSubstrateNode :
    public std::unary_function<xmlpp::Node*, std::map<int, int>::value_type>
  {
    tauSpeciesCatalog& rSpeciesCatalog;
    
  public:
    parseReactionSubstrateNode(tauSpeciesCatalog& rTauSpeciesCatalog) :
      rSpeciesCatalog(rTauSpeciesCatalog)
    {}

    result_type
    operator()(xmlpp::Node* pReactionSubstrateNode) throw(std::exception)
    {
      xmlpp::Element* pReactionSubstrateElt
	= utl::dom::mustBeElementPtr(pReactionSubstrateNode);

      std::string substrateName
	= utl::dom::mustGetAttrString(pReactionSubstrateElt,
				      tauEltName::reactionSubstrate_nameAttr);
      int speciesNdx
	= rSpeciesCatalog.mustFindSpecies(pReactionSubstrateElt,
					  substrateName);
      int multiplicity
	= utl::dom::mustGetAttrInt(pReactionSubstrateElt,
				   tauEltName::reactionSubstrate_multAttr);

      return std::make_pair(speciesNdx,
			    multiplicity);
    }
  };

  // Parse a reaction product species and its multiplicity into the reaction
  // deltas, a column in the stoichiometric matrix.  The multiplicities of all
  // the reactant species have already been inserted as negtive deltas.
  class parseReactionProductNode :
    public std::unary_function<xmlpp::Node*, void>
  {
    tauSpeciesCatalog& rSpeciesCatalog;
    parserReaction& rReaction;

  public:
    parseReactionProductNode(tauSpeciesCatalog& rTauSpeciesCatalog,
			     parserReaction& rParserReaction) :
      rSpeciesCatalog(rTauSpeciesCatalog),
      rReaction(rParserReaction)
    {}

    void
    operator()(xmlpp::Node* pReactionProductNode) throw(std::exception)
    {
      xmlpp::Element* pReactionProductElt
	= utl::dom::mustBeElementPtr(pReactionProductNode);

      std::string productSpeciesName
	= utl::dom::mustGetAttrString(pReactionProductElt,
				      tauEltName::reactionProduct_nameAttr);
      int speciesNdx
	= rSpeciesCatalog.mustFindSpecies(pReactionProductElt,
					  productSpeciesName);
      int multiplicity
	= utl::dom::mustGetAttrInt(pReactionProductElt,
				   tauEltName::reactionProduct_multAttr);

      // See if the species is not a substrate by attempting to insert
      // the species delta pair.  This will give the desired
      // result if successful.
      std::pair<std::map<int,int>::iterator, bool> insertResult
	= rReaction.deltas.insert(std::make_pair(speciesNdx,
						 multiplicity));

      // If the species is a substrate, then its current delta
      // should be the negative of its multiplicity as a reactant.
      // Hence, we just add its multiplicity as a product.
      if(! insertResult.second)
	{
	  insertResult.first->second += multiplicity;
	}
    }
  };

  class initDeltaFromSubstrate :
    public std::unary_function<std::map<int,int>::value_type,
    std::map<int,int>::value_type>
  {
  public:
    result_type
    operator()(const argument_type& rSubstrateMultPair) const
    {
      return std::make_pair(rSubstrateMultPair.first,
			    - rSubstrateMultPair.second);
    }
  };

  class parseReaction :
    public std::unary_function<xmlpp::Node*, parserReaction>
  {
    tauSpeciesCatalog& rSpeciesCatalog;
  public:
    parseReaction(tauSpeciesCatalog& rTauSpeciesCatalog) :
      rSpeciesCatalog(rTauSpeciesCatalog)
    {}
  
    parserReaction
    operator()(xmlpp::Node* pReactionNode) throw(std::exception)
    {
      parserReaction reaction;

      xmlpp::Element* pReactionElt
	= utl::dom::mustBeElementPtr(pReactionNode);

      // Parse the substrates and their multiplicities.
      xmlpp::Node::NodeList substrateNodes
	= pReactionElt->get_children(tauEltName::reactionSubstrate);
      std::transform(substrateNodes.begin(),
		     substrateNodes.end(),
		     std::inserter(reaction.substrates,
				   reaction.substrates.begin()),
		     parseReactionSubstrateNode(rSpeciesCatalog));

      // This is a patch to cover up the change in format for dump of
      // reactions in moleculizer state dumps.  That change was made
      // to support SBML, I think.
      //
      // Initialize the deltas with the negations of the substrate
      // multiplicities.  (Now, every time I see "substrate" it sets
      // off alarm bells as per Gillespie.)
      std::transform(reaction.substrates.begin(),
		     reaction.substrates.end(),
		     std::inserter(reaction.deltas,
				   reaction.deltas.begin()),
		     initDeltaFromSubstrate());

      // Now, complete the patch by parsing product nodes by adding their
      // multiplicities to the negative deltas coming from reactants.
      xmlpp::Node::NodeList productNodes
	= pReactionElt->get_children(tauEltName::reactionProduct);
      std::for_each(productNodes.begin(),
		    productNodes.end(),
		    parseReactionProductNode(rSpeciesCatalog,
					     reaction));

      // Parse the reaction rate.
      xmlpp::Element* pRateElt
	= utl::dom::mustGetUniqueChild(pReactionElt,
				       tauEltName::rate);
      reaction.rate
	= utl::dom::mustGetAttrDouble(pRateElt,
				      tauEltName::rate_valueAttr);

      return reaction;
    }
  };

  class lookUpDumpableSpecies :
    public std::unary_function<xmlpp::Node*, int>
  {
    tauSpeciesCatalog& rSpeciesCatalog;
  public:
    lookUpDumpableSpecies(tauSpeciesCatalog& rTauSpeciesCatalog) :
      rSpeciesCatalog(rTauSpeciesCatalog)
    {}

    int
    operator()(xmlpp::Node* pSpeciesRefNode) throw(std::exception)
    {
      xmlpp::Element* pSpeciesRefElt
	= utl::dom::mustBeElementPtr(pSpeciesRefNode);

      std::string speciesName
	= utl::dom::mustGetAttrString(pSpeciesRefElt,
				      tauEltName::speciesRef_nameAttr);

      return rSpeciesCatalog.mustFindSpecies(pSpeciesRefElt,
					     speciesName);
    }
  };

  // Parses a "dumpable" node, which contains the name of the dumpable,
  // destined to be a column header in a dump file, along with a list
  // of one or more species.
  class parseDumpableNode :
    public std::unary_function<xmlpp::Node*, tauDumpableCatalog::value_type>
  {
    tauSpeciesCatalog& rSpeciesCatalog;
  public:
    parseDumpableNode(tauSpeciesCatalog& rTauSpeciesCatalog) :
      rSpeciesCatalog(rTauSpeciesCatalog)
    {}

    result_type
    operator()(xmlpp::Node* pDumpableNode) throw(std::exception)
    {
      tauDumpable dumpable;
    
      xmlpp::Element* pDumpableElt
	= utl::dom::mustBeElementPtr(pDumpableNode);

      dumpable.name
	= utl::dom::mustGetAttrString(pDumpableElt,
				      tauEltName::dumpable_nameAttr);

      xmlpp::Node::NodeList speciesRefNodes
	= pDumpableElt->get_children(tauEltName::speciesRef);
      std::transform(speciesRefNodes.begin(),
		     speciesRefNodes.end(),
		     back_inserter(dumpable),
		     lookUpDumpableSpecies(rSpeciesCatalog));

      return std::make_pair(dumpable.name,
			    dumpable);
    }
  };

  class parseDumpableRefNode :
    public std::unary_function<xmlpp::Node*, tauDumpable>
  {
    const tauDumpableCatalog& rDumpableCatalog;
  public:
    parseDumpableRefNode(const tauDumpableCatalog& rTauDumpableCatalog) :
      rDumpableCatalog(rTauDumpableCatalog)
    {}
  
    const tauDumpable&
    operator()(xmlpp::Node* pDumpableRefNode) throw(std::exception)
    {
      xmlpp::Element* pDumpableRefElt
	= utl::dom::mustBeElementPtr(pDumpableRefNode);

      std::string dumpableName
	= utl::dom::mustGetAttrString(pDumpableRefElt,
				      tauEltName::dumpableRef_nameAttr);

      return rDumpableCatalog.mustFindDumpable(pDumpableRefElt,
					       dumpableName);
    }
  };

  class parseDumpStreamNode :
    public std::unary_function<xmlpp::Node*, dumpStream>
  {
    const tauDumpableCatalog& rDumpableCatalog;
  public:
    parseDumpStreamNode(const tauDumpableCatalog& rTauDumpableCatalog) :
      rDumpableCatalog(rTauDumpableCatalog)
    {}

    dumpStream
    operator()(xmlpp::Node* pDumpStreamNode) throw(std::exception)
    {
      dumpStream result;
    
      xmlpp::Element* pDumpStreamElt
	= utl::dom::mustBeElementPtr(pDumpStreamNode);

      // Parse the file name, which seems to be given by the "name"
      // attribute of the dumpStream.
      std::string fileName
	= utl::dom::mustGetAttrString(pDumpStreamElt,
				      tauEltName::dumpStream_nameAttr);

      // Parse the dumpables that are supposed to appear as columns
      // in the output file.
      xmlpp::Node::NodeList dumpableRefNodes
	= pDumpStreamElt->get_children(tauEltName::dumpableRef);

      std::transform(dumpableRefNodes.begin(),
		     dumpableRefNodes.end(),
		     std::back_inserter(result),
		     parseDumpableRefNode(rDumpableCatalog));

      result.mustSetStream(pDumpStreamElt,
			   fileName);
      return result;
    }
  };

  tauApp::tauApp(int argc,
		 char** argv,
		 xmlpp::Document* pDoc) throw(std::exception) :
    now(0.0),
    dumpTime(0.0),
    dumpInterval(1.0),
    stopTime(1.0),
    epsilon(1.0),
    seedString("A seedy string."),
    molarConst(1.0),
    pAdaptor(0)
  {
    xmlpp::Element* pRootElt
      = pDoc->get_root_node();

    xmlpp::Element* pSpeciesElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::species);

    // Parse the species names and populations.
    tauSpeciesCatalog speciesCatalog;
    xmlpp::Node::NodeList speciesEntryNodes
      = pSpeciesElt->get_children(tauEltName::speciesEntry);

    std::for_each(speciesEntryNodes.begin(),
		  speciesEntryNodes.end(),
		  parseSpeciesEntry(speciesCatalog,
				    populations));

    xmlpp::Element* pReactionsElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::reactions);

    // Parse the reactions.
    //
    // For now, I'm sticking here with the substrates/deltas format, though
    // I will probably need to switch to substrates/products format.
    //
    // As noted above, I could probably do away with parserReaction by making
    // several passes over the reactions in the doc?
    xmlpp::Node::NodeList reactionNodes
      = pReactionsElt->get_children(tauEltName::reaction);

    std::vector<parserReaction> reactions;
    std::transform(reactionNodes.begin(),
		   reactionNodes.end(),
		   std::back_inserter(reactions),
		   parseReaction(speciesCatalog));

    xmlpp::Element* pDumpablesElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::dumpables);

    // Parse the dumpables.  These are used to cook up dumpStreams, which
    // are all that the program knows.
    xmlpp::Node::NodeList dumpableNodes
      = pDumpablesElt->get_children(tauEltName::dumpable);
    tauDumpableCatalog dumpables;
    std::transform(dumpableNodes.begin(),
		   dumpableNodes.end(),
		   std::inserter(dumpables,
				 dumpables.begin()),
		   parseDumpableNode(speciesCatalog));

    // Parse the dump-streams.
    xmlpp::Element* pDumpStreamsElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::dumpStreams);

    xmlpp::Node::NodeList dumpStreamNodes
      = pDumpStreamsElt->get_children(tauEltName::dumpStream);

    // The name of a dumpStream is really a file name.  In the old parser code,
    // the file was opened when the name of the stream was read.
    std::transform(dumpStreamNodes.begin(),
		   dumpStreamNodes.end(),
		   std::back_inserter(dumpStreams),
		   parseDumpStreamNode(dumpables));

    // Parse the various other constants.
    //
    // Volume, which is stored in the same way as moleculizer.
    xmlpp::Element* pVolumeElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::volume);
    double volume
      = utl::dom::mustGetAttrDouble(pVolumeElt,
				    tauEltName::volume_litersAttr);
    molarConst = volume * fnd::avogadrosNumber;

    // Stop-time
    xmlpp::Element* pStopTimeElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::stopTime);
    stopTime
      = utl::dom::mustGetAttrDouble(pStopTimeElt,
				    tauEltName::stopTime_secondsAttr);

    // Dump interval.
    xmlpp::Element* pDumpIntervalElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::dumpInterval);
    dumpInterval
      = utl::dom::mustGetAttrDouble(pDumpIntervalElt,
				    tauEltName::dumpInterval_secondsAttr);

    // Epsilon
    xmlpp::Element* pEpsilonElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::epsilon);
    epsilon
      = utl::dom::mustGetAttrDouble(pEpsilonElt,
				    tauEltName::epsilon_valueAttr);

    pAdaptor = new ckStepAdaptor(epsilon);

    // Seed string.  Hashed into an integer by the main program?
    xmlpp::Element* pSeedElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     tauEltName::seed);
    seedString
      = utl::dom::mustGetAttrString(pSeedElt,
				    tauEltName::seed_valueAttr);

    // These two steps use the volume as well as the parsed reactions,
    // so they have to come at the end here.
    
    // Generate the "stoichiometric matrix", which is actually
    // coded as a linear polynomial function.
    makeCountConverter(reactions);

    // Generate the polynomials that give the derivatives of the reaction
    // amounts in terms of their values, the underlying (transformed)
    // differential equation.
    makeDerivatives(reactions);

  }
}

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

#include <algorithm>
#include "plex/recognizer.hh"
#include "plex/plexFamily.hh"
#include "plex/plexMap.hh"
#include "plex/plexQuery.hh"
#include "mzr/util.hh"

namespace plx
{
  recognizer::~recognizer(void)
  {
    // The deleter for map values should work with multimaps, too.
    for_each(plexHasher.begin(),
	     plexHasher.end(),
	     mzr::deleteMapValue<int, plexFamily*>());
  }

  // Function class for finding a plex species among those with the same
  // hash value.  This version tracks the isomorphism with the
  // plexFamily's paradigm.
  class thisPlexFamilyIso : std::unary_function<std::pair<const int, plexFamily*>, bool>
  {
    const plex& rPlex;
    plexIsoPair& rIsomorphism;

  public:
    thisPlexFamilyIso(const plex& rPlexToClassify,
		      plexIsoPair& rIso) :
      rPlex(rPlexToClassify),
      rIsomorphism(rIso)
    {}

    bool operator()(const std::pair<const int, plexFamily*>& rHashEntry)
    {
      return reportIsoSearch(rPlex,
			     rHashEntry.second->getParadigm(),
			     rIsomorphism).findIso();
    }
  };

  // Function class for finding a plex species among those with the same
  // hash value.  This one doesn't do tracking of the isomorphism.
  class isThisPlexFamily : std::unary_function<std::pair<const int, plexFamily*>, bool>
  {
    const plex& rPlex;
  public:
    isThisPlexFamily(const plex& rPlexToClassify) :
      rPlex(rPlexToClassify)
    {}

    bool operator()(const std::pair<const int, plexFamily*>& rHashEntry) const
    {
      return plexIsoSearch(rPlex, rHashEntry.second->getParadigm()).findIso();
    }
  };

  bool
  recognizer::unify(const plex& aPlex,
		    plexFamily*& rpFamily,
		    plexIsoPair* pIso)
  {
    // Check the cache to see if this plex has ever been recognized
    // before.
    std::pair<std::map<plex, recognition>::iterator, bool> insertResult
      = recognizedCache.insert(std::make_pair(aPlex,
					      recognition()));

    // We use these references to correct the entry in the
    // cache if the insert (of the invalid recognition)
    // succeeds.
    plexFamily*& rFamilyPtr = insertResult.first->second.pPlexFamily;
    plexIsoPair& rIso = insertResult.first->second.iso;

    // If the insert succeeded, then the plex is either totally
    // unknown or is isomorphic to a known plex.
    bool familyIsNew = false;
    if(insertResult.second)
      {
	// Get the plex's hash value.
	int plexHashValue = aPlex.hashValue();

	// Get iterator range to all plex iso classes with this hash
	// value.
	std::pair<std::multimap<int, plexFamily*>::iterator,
	  std::multimap<int, plexFamily*>::iterator> rangeIterPair
	  = plexHasher.equal_range(plexHashValue);

	// Scan the plex isomorphism classes having the same hash value
	// as this plex for the correct isomorphism class of this plex.
	std::multimap<int, plexFamily*>::iterator iEntry
	  = find_if(rangeIterPair.first,
		    rangeIterPair.second,
		    thisPlexFamilyIso(aPlex,
				      rIso));

	if(iEntry == rangeIterPair.second)
	  {
	    // We have never seen the plex before, so we have to construct
	    // its plexFamily.
	    familyIsNew = true;
	  
	    // Construct a new plexFamily.
	    rFamilyPtr = new plexFamily(aPlex,
					rPlexUnit);

	    // The given plex is the paradigm of the new plexFamily,
	    // so the isomorphism is the identity.
	    rIso = plexIsoPair::makeIdentity(aPlex.mols.size(),
					     aPlex.bindings.size());

	    // Rememember this family, in case we ever see it again.
	    plexHasher.insert(std::make_pair(plexHashValue,
					     rFamilyPtr));
	  }
	else
	  {
	    // The plex belongs to a family that we've already seen before,
	    // but was not in the cache.
	    rFamilyPtr = iEntry->second;
	  }
      }

    // Return plexFamily pointer that was installed in the map.
    rpFamily = rFamilyPtr;

    // Optionally return (i.e. copy out) the isomorphism from the
    // plex to the family's paradigm plex.
    if(pIso) *pIso = rIso;

    // Return whether plexFamily was recognized for the first time.
    return familyIsNew;
  }

  // Finds the isomorphism class (plexFamily) of a given plex.  Does not
  // track the isomorphism with the plexFamily's paradigm.
  //
  // Note that this routine doesn't let you know if the plex was new or
  // old.  That feature would be easy to add, and it would make it
  // possible to notify the user if (s)he'd defined the same plex under
  // two different names.
  //
  // This is the most commonly used form of recognition. It is used at
  // "run time," after all definitions and specifications are complete.
  // It recognizes omniPlexes in the recognized plexFamily.
  plexFamily*
  recognizer::operator()(const plex& aPlex)
  {
    // Do the basic recognition, without initialization of the plexFamily.
    plexFamily* pPlexFamily;
    if(unify(aPlex,
	     pPlexFamily,
	     0))
      {
	// Connect the plexFamily to all its features.
	pPlexFamily->connectToFeatures();
      }

    return pPlexFamily;
  }

  // Does recognition, including finding active subcomplexes (omniPlexes),
  // and reports the isomorphism with the plexFamily's paradigm.  This
  // version is for use at run time, and assumes that all omniplexes are
  // already specified.
  //
  // This is used in the decomposition reaction.
  plexFamily*
  recognizer::operator()(const plex& aPlex,
			 plexIsoPair& rIso)
  {
    // Do the basic recognition, without initialization of the plexFamily.
    plexFamily* pPlexFamily;
    if(unify(aPlex,
	     pPlexFamily,
	     &rIso))
      {
	// Connect the plexFamily to all its features.
	pPlexFamily->connectToFeatures();
      }

    return pPlexFamily;
  }

  class insertFamilySpecies :
    public std::unary_function<std::map<int, plexFamily*>::value_type, void>
  {
    xmlpp::Element* pExplicitSpeciesElt;
    double molFact;
  public:
    insertFamilySpecies(xmlpp::Element* pExplicitSpeciesElement,
			double molarFactor) :
      pExplicitSpeciesElt(pExplicitSpeciesElement),
      molFact(molarFactor)
    {}

    void
    operator()(const argument_type& rHasherEntry) const throw(std::exception)
    {
      const plexFamily* pFamily = rHasherEntry.second;
      pFamily->insertSpecies(pExplicitSpeciesElt,
			     molFact);
    }
  };

  void
  recognizer::insertSpecies(xmlpp::Element* pExplicitSpeciesElt,
			    double molarFactor) const
    throw(std::exception)
  {
    std::for_each(plexHasher.begin(),
		  plexHasher.end(),
		  insertFamilySpecies(pExplicitSpeciesElt,
				      molarFactor));
  }
}

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

#ifndef MOL_H
#define MOL_H

#include <vector>
#include <string>
#include <sstream>
#include "mzr/species.hh"
#include "mzr/feature.hh"
#include "plex/prm.hh"
#include "plex/plexSpecies.hh"
#include "mol/molXcpt.hh"
#include "mol/bindingSite.hh"

namespace bnd
{
  class mol :
    public mzr::feature<plx::plexSpecies, plx::plexMolSpec>
  {
    // Note that there is no public way to add a binding site after
    // construction.  This is because some descendant classes call
    // getDefaultSiteParams during their construction, so that all
    // the sites have to be in place by the end of this base-class
    // constructor.
    class addSite :
      public std::unary_function<const bindingSite&, void>
    {
      mol& rMol;
    public:
      addSite(mol& rTargetMol) :
	rMol(rTargetMol)
      {}

      void
      operator()(const bindingSite& rSite) const
      {
	int siteNdx = rMol.sites.size();
	if(! rMol.siteNameToNdx.insert
	   (std::pair<std::string, int>(rSite.getName(),
					siteNdx)).second)
	  throw duplicateSiteNameXcpt(rMol,
				      rSite.getName());

	rMol.sites.push_back(rSite);
      }
    };
    friend class addSite;
      
  protected:
  
    std::string name;

    // Binding sites are accessed and referred to
    // by index, except by the user.
    std::vector<bindingSite> sites;

    // This is for user/command use.
    std::map<std::string, int> siteNameToNdx;

  public:

    // Note that there is no public way to add a binding site after
    // construction.  This is because some descendant classes call
    // getDefaultSiteParams during their construction, so that all
    // the sites have to be in place by the end of this base-class
    // constructor.
    mol(const std::string& rName,
	const std::vector<bindingSite>& rSites) :
      name(rName)
    {
      for_each(rSites.begin(),
	       rSites.end(),
	       addSite(*this));
    }

    virtual
    ~mol(void)
    {}
  
    const std::string&
    getName(void) const
    {
      return name;
    }

    int
    getSiteCount(void) const
    {
      return sites.size();
    }

    // I don't think this can be const; or at least, there has to
    // be a non-const version of it.
    //
    // Added later: this function couldd go away, now that bindingSite
    // inherits from feature.
    mzr::feature<plx::plexSpecies, plx::plexSiteSpec>&
    getSiteFeature(int siteNdx)
    {
      return sites[siteNdx];
    }

    // Convert site name to site index.  Badly named, given the new
    // status of binding sites as real entities.
    bool
    findSite(const std::string& rName,
	     int& rSiteNdx) const
    {
      std::map<std::string, int>::const_iterator iEntry
	= siteNameToNdx.find(rName);

      bool siteNameFound
	= (iEntry != siteNameToNdx.end());

      if(siteNameFound) rSiteNdx = iEntry->second;

      return siteNameFound;
    }

    int
    mustGetSiteNdxByName(const xmlpp::Node* pRequestingNode,
			 const std::string& rName) const
      throw(unknownSiteXcpt)
    {
      int siteNdx = -1;
      if(! findSite(rName,
		    siteNdx))
	throw unknownSiteXcpt(pRequestingNode,
			      rName);
      return siteNdx;
    }

    // For adding shapes to binding site.  Returns 0 if no such name.
    // This is temporary, in support of the "site" scripting command.
    bindingSite*
    getSiteByName(const std::string& rName)
    {
      std::map<std::string, int>::iterator iEntry
	= siteNameToNdx.find(rName);

      if(siteNameToNdx.end() == iEntry) return 0;
      else return &(sites[iEntry->second]);
    }

    const bindingSite*
    getSiteByName(const std::string& rName) const
    {
      std::map<std::string, int>::const_iterator iEntry
	= siteNameToNdx.find(rName);

      if(siteNameToNdx.end() == iEntry) return 0;
      else return &(sites[iEntry->second]);
    }

    // There are some situations where we need to get a site
    // by index, such as from a plexSiteSpec, such as in the
    // plex command alloSiteCmd.  Given this, it might be easier/
    // better just to expose the vector of binding sites.
    bindingSite&
    getSiteByNdx(int siteNdx)
    {
      return sites[siteNdx];
    }

    const bindingSite&
    getSiteByNdx(int siteNdx) const
    {
      return sites[siteNdx];
    }

    // This has to be non-const for various reasons:
    // 
    // During mol definition, to set the allosteric shapes of the sites
    // when defining an allosteric state.
    //
    // If the param is encountered for the first time (as after a new
    // modification to the mol) this routine "interns" the new state.
    virtual std::vector<siteParam>&
    allostery(molParam param) = 0;

    // This is support for the construction of the default state of
    // a plex, which starts from the vector of default states of its
    // mols.
    virtual molParam
    getDefaultParam(void) const = 0;

    // Returns vector of default site shape pointers for the sites.
    std::vector<siteParam>
    getDefaultSiteParams(void) const
    {
      std::vector<siteParam> defaultParams(getSiteCount());

      transform(sites.begin(),
		sites.end(),
		defaultParams.begin(),
		mem_fun_ref(&bindingSite::getDefaultParam));

      return defaultParams;
    }

    // This is how mols in species of complexes will describe their state.
    // Mols that have only one state do nothing and return null pointer.
    virtual xmlpp::Element*
    insertInstanceState(xmlpp::Element* pInstanceStatesElt,
			int molInstanceNdx,
			molParam param) const
    {
      return 0;
    }

    // This is how an instance name for a mol in a complex is automatically
    // generated from its index in the paradigm plex.  This virtual function
    // arranges that different kinds of mols can distinguish themselves
    // slightly, for the user's convenience.
    //
    // This routine is used by plexes to give their structure and by
    // the this mol to convert a state pointer and and instance index
    // into an instance state description (a "mod-map" for a mod-mol).
    //
    // This has to be done for all automatically generated plex species, and for
    // the time being, I'm also going to do it for user-specified plex species.
    // That is, when state is dumped user-specified plex species will be dumped
    // with these automatically generated instance names, rather than any
    // instance names the user may have given.  (Actually, a user could use
    // different instance names when specifying the plex in different contexts.)
    virtual std::string
    genInstanceName(int molInstanceNdx) const
    {
      std::ostringstream oss;
      oss << "mol_"
	  << molInstanceNdx;
      return oss.str();
    }

    // This is how a mol inserts itself into the state document (assuming now
    // that the state document is going to generate itself from scratch each
    // time it's called for.)
    virtual xmlpp::Element*
    insertElt(xmlpp::Element* pMolsElt) const throw(std::exception) = 0;

    // Note that the target vector here must be of the correct size,
    // which can be obtained with getSiteCount().
    void
    getSiteNames(std::vector<std::string>& rSiteNames) const
    {
      transform(sites.begin(),
		sites.end(),
		rSiteNames.begin(),
		mem_fun_ref(&bindingSite::getName));
    }

    std::vector<std::string>
    getSiteNames(void) const
    {
      std::vector<std::string> siteNames(getSiteCount());
      getSiteNames(siteNames);
      return siteNames;
    }
  };

  // Base class for 
  template<class stateClass>
  class stateMol :
    public mol
  {
  protected:
    // This pointer is null in this base class.
    //
    // Usually, the defaultState is somewhere to be found in the
    // stateMols's member data, varying from class to class.
    const stateClass* pDefaultState;

    // Having this vector pre-extracted from the binding sites accelerates
    // interning new states in allomols and the allostery function in
    // simple mols.
    std::vector<siteParam> defaultSiteParams;
  public:
    typedef stateClass molStateClass;
  
    stateMol(const std::string& rName,
	     const std::vector<bindingSite>& rSites) :
      mol(rName,
	  rSites),
      pDefaultState(0)
    {
      defaultSiteParams = getDefaultSiteParams();
    }

    // This is used to get the default state as a starting point in the
    // construction of new states in various situations around the
    // program.
    //
    // Note that it returns a pointer to the particular class of state
    // for this mol, while the more generic getDefaultParam returns a
    // pointer to the base class molState (see plex/prm.h).
    //
    // It is badly named.
    const stateClass*
    getDefaultState(void) const
    {
      return pDefaultState;
    }

    virtual std::vector<siteParam>&
    allostery(molParam param)
    {
      return defaultSiteParams;
    }

    virtual molParam
    getDefaultParam(void) const
    {
      return getDefaultState();
    }
  };

  /*! \ingroup plexSpeciesGroup

  \brief Template base class for mols with many states, exhibiting
  allostery. */
  template<class stateClass>
  class alloMol :
    public stateMol<stateClass>
  {
  protected:

    // This mapping is where mol states are interned.
    typedef std::map<stateClass, std::vector<siteParam> > alloMapType;
    alloMapType alloMap;

  public:
    // The default state can't be interned (to set
    // pInternedDefaultState) until all the sites have been added to the
    // mol, so that the vector siteParams associated to the default
    // state can be obtained.
    alloMol(const std::string& rName,
	    const std::vector<bindingSite>& rSites,
	    const stateClass& rDefaultState) :
      stateMol<stateClass>(rName,
			   rSites)
    {
      pDefaultState = internState(rDefaultState);
    }

    // Sometimes, the default state can't be constructed until the rest of the
    // mol is finished. (In modMol, the map from modification site names
    // to modification site indices is needed before the modifications
    // can be put into a vector.)  This constructor initializes pDefaultState
    // to 0, leaving it for later.
    alloMol(const std::string& rName,
	    const std::vector<bindingSite>& rSites) :
      stateMol<stateClass>(rName,
			   rSites)
    {}

    // phosMol, for example, overrides this to check that the state
    // object has the right number of phosphorylation sites.
    virtual const stateClass*
    internState(const stateClass& rState)
    {
      typename alloMapType::iterator iStateEntry
	= alloMap.find(rState);

      if(alloMap.end() == iStateEntry)
	{
	  iStateEntry
	    = alloMap.insert(std::make_pair(rState,
					    defaultSiteParams)).first;
	}

      return &(iStateEntry->first);
    }

    // This could be static, except that we want to report the
    // name of the mol (at least) in the exception.
    const stateClass&
    externState(molParam prm) const
    {
      const stateClass* pOurState = dynamic_cast<const stateClass*>(prm);

      if(! pOurState) throw badMolParamClassXcpt(getName());

      return *pOurState;
    }

    std::vector<siteParam>&
    allostery(molParam param)
    {
      const stateClass& rState = externState(param);

      // Look up the state in the allosteryMap, installing the
      // defaultSiteParams if not found.
      typename alloMapType::iterator iStateEntry
	= alloMap.find(rState);

      if(alloMap.end() == iStateEntry)
	{
	  iStateEntry
	    = alloMap.insert(std::make_pair(rState,
					    defaultSiteParams)).first;
	}

      return iStateEntry->second;
    }
  };
}

#endif // MOL_H

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

#ifndef MODMIXIN_H
#define MODMIXIN_H

// These classes mix in the ability to add modifications at named sites to
// mols.  The mol base class is a fundamental class for dealing with
// essentially symmetric binding between complexes.  Modifications are
// for dealing with binding with entities that are treated in the
// simulation as featureless, except for molecular weight, typically
// small molecules.  Binding to ATP or phosphate are important cases
// in our current simulations.
//
// There is one mixin for the mol class and another for the molState
// class.

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "mzr/mzrXcpt.hh"
#include "mzr/util.hh"
#include "sampleDist/sampleDist.hh"
#include "mol/modification.hh"
#include "mol/molDomParse.hh"

namespace bnd
{
  class modMol;

  // This is the contribution of this mixin to the overall mol state.
  // 
  // It is also used as the "core" of a query in which modifications
  // match themselves and the null pointer represents a wildcard.
  class modStateMixin :
    public std::vector<const modification*>
  {
  public:
    modStateMixin(const modification* pMod,
		  int count) :
      std::vector<const modification*>(count, pMod)
    {}
  
    modStateMixin(const std::vector<const modification*>& rModifications) :
      std::vector<const modification*>(rModifications)
    {}

    double
    totalWeightDelta(void) const;

    // Gives the number of modification sites with the given modification.
    int
    modCount(const modification* pModToCount) const;

    // Same as above, but only counts sites that are masked true.
    int
    modCount(const modification* pModToCount,
	     const std::vector<bool>& rMask) const;

    // Changes a randomly selected site with modification pFromMod to
    // have modification pToMod.  Expect this to be used, e.g. to
    // phosphorylate an unphosphorylated site.
    //
    // Note that this could just as well be used to dephosphorylate
    // a phosphorylated site.
    // 
    // Returns the index of the modified site or -1 if no site had
    // modification fromMod.
    //
    // Do I need these any more?
    //   int
    //   pushMod(const modification* pFromMod,
    // 	  const modification* pToMod);

    //   // Same as above, but acts only on sites masked true.
    //   int
    //   pushMod(const modification* pFromMod,
    // 	  const modification* pToMod,
    // 	  std::vector<bool>& rMask);
  };

  // This mixin adds the ability to access modification states.
  //
  // An elaboration here would be for each modification site to have a list of
  // acceptable modifications.
  class modMolMixin
  {
  public:
    // These mappings allow the modification sites of the mol to be accessed
    // by name or by index.
    std::map<std::string, int> modSiteNameToNdx;
    std::vector<std::string> modSiteNames;

    int
    getModSiteNdx(const std::string& rModSiteName) const;

    // This should be an method of modMol, obviating the stupid pModMol
    // argument.
    int
    mustGetModSiteNdx(xmlpp::Node* pRequestingNode,
		      const std::string& rModSiteName,
		      modMol* pModMol) const
      throw(unknownModSiteXcpt);

    int
    modSiteCount(void) const
    {
      return modSiteNames.size();
    }

    // I will use a map (encoded into XML) from modification site names
    // to modification names in the following circumstances:
    //
    // 1. To simultaneously define the modification sites and give their
    // default modifications.
    //
    // 2. To specify a single modification state by giving its differences
    // from the default state.  That is, any blank entries are treated
    // as representing the default state.
    //
    // 3. To specify a collection of modification states by specifying
    // the modifications at some of the sites.  Blank entries are
    // treated as wildcards; i.e. no filtering is done on unspecified
    // modification sites.  I will have to be able to "or" these
    // together.

    // Constructs the mol's map from modification site names to modification
    // site indices.  This accomplishes part of goal 1 above; the rest
    // of goal 1 comes from also using indexModMap below.
    modMolMixin
    (const std::map<std::string, const modification*>& rDefaultModMap);


    // Constructs a modStateMixin by doing substitutions on a provided
    // modStateMixin.  This would be used to construct a state from
    // a map and the default state.  This accomplishes goal 2 above.
    modStateMixin
    substituteModMap(const std::map<std::string, const modification*>& rModMap,
		     const modStateMixin& rSourceStateMixin);

    // Uses the map from modification site names to modification
    // site indices to convert a map from modification site names
    // to modifications into a vector of modification*'s.
    //
    // If any modification site is unmapped, then the corresponding
    // pointer in the modStateMixin will be null.  This makes this
    // function useful for constructing "regexp" modStateMixin's for use
    // with modStateMatch below.  This is the bulk of goal 3.
    //
    // Another use is in the construction of the modStateMixin
    // part of the default state.
    modStateMixin
    indexModMap(const std::map<std::string, const modification*>& rModMap);

    // This predicate can be used to test modStateMixin's against a
    // "regexp" modStateMixin.  If a null pointer occurs as an entry in
    // the "regexp," then it matches any modification.  Otherwise the
    // modifications (modification*'s) must match exactly.
    //
    // The function indexModMap above should be useful for constructing
    // these "regexps", since null pointers will naturally appear at
    // unmapped indices.  Together with indexModMap above, this
    // accomplishes goal 3.
    class modStateMatch :
      public std::unary_function<modStateMixin, bool>
    {
      const modStateMixin& rMatch;
    public:
      modStateMatch(const modStateMixin& rMatchStateMixin) :
	rMatch(rMatchStateMixin)
      {}

      bool
      operator()(const modStateMixin& rMixinToTest) const;
    };
  };
}

#endif // MODMIXIN_H

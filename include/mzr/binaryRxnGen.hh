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

#ifndef BINARYRXNGEN_H
#define BINARYRXNGEN_H

/*! \file
  \ingroup rxnGenGroup
  \brief Defines reaction generator template for "binary" reaction families.

  Three templates are given in all.  The first two support the last,
  which is the only one intended for direct instantiation.  This last
  template is useful as a base class for reaction generators for any
  reactionFamily that is sensitive to the creation of new species
  reported by two members of feature.

  For example, the dimerizeFam family of dimerization reactions is
  sensitive to two free binding site features.  Whenever a complex
  species that displays one of these free binding sites is created,
  the dimerizeFam must create new dimerizations between the new
  species of complex and all the species displaying a complementary
  site.  */

#include <functional>
#include "mzr/featureContext.hh"

namespace mzr
{
  // These templates implement reaction generation for binary reaction
  // families.  The family is binary in the sense that it can pay attention
  // to two features and generate all the appropriate new reactions when
  // either feature reports a new species.
  //
  // This is done by implementing two reaction generators, one for each
  // of the two features to which the reaction family is connected.
  //
  // The two generators are the same, but are used in complementary
  // fashion in the binary family.  Below, binaryRxnGenPair
  // implements a complementary pair for use in creating new families of
  // binary reactions.

  /*! \ingroup rxnGenGroup
    \brief Support template for the binaryRxnGenPair template.

    The full reaction generator for a binary reaction has two mirror-image
    "sub-generators."  Each sees to adding reactions coming from one of
    the features to which the reaction family is sensitive. */
  template<class newContextClass, class otherFeatureClass>
  class binaryRxnGen :
    public newContextClass::rxnGen
  {

    typedef typename otherFeatureClass::context otherContextClass;

    otherFeatureClass& rOtherFtr;

    class forOneOtherContext :
      public std::unary_function<otherContextClass, void>
    {
      const newContextClass& rNew;
      const binaryRxnGen& rRxnGen;
    public:
      forOneOtherContext(const newContextClass& rNewContext,
			 const binaryRxnGen& rThisRxnGen) :
	rNew(rNewContext),
	rRxnGen(rThisRxnGen)
      {}

      // Make and install the parameters resulting from combining the new
      // context rNew with one context coming from the other rxnGen.
      void operator()(const otherContextClass& rOtherContext) const
      {
	rRxnGen.makeBinaryReactions(rNew,
				    rOtherContext);
      }
    };
    friend class forOneOtherContext;

  public:
    binaryRxnGen(otherFeatureClass& rOtherFeature) :
      rOtherFtr(rOtherFeature)
    {}

    void
    makeReactions(const newContextClass& rContext) const
    {
      for_each(rOtherFtr.contexts.begin(),
	       rOtherFtr.contexts.end(),
	       forOneOtherContext(rContext,
				  *this));
    }
    
    virtual void
    makeBinaryReactions(const newContextClass& rNewContext,
			const otherContextClass& rOtherContext) const = 0;
  
  };

  /*! \ingroup rxnGenGroup
    \brief Template reaction generator for binary reaction.

    Basically, two complementary reaction generators that act in tandem
    are implemented, one for each of the features.  */
  template<class leftFeatureClass, class rightFeatureClass>
  class binaryRxnGenPair
  {
    typedef typename leftFeatureClass::context leftContextClass;
    typedef typename rightFeatureClass::context rightContextClass;

    class leftRxnGenClass :
      public binaryRxnGen<leftContextClass, rightFeatureClass>
    {
      // Normally, I would use a reference here, but that prevents
      // this class from being acceptable as a container element.
      // Some containers use assignment, instead of copy constructors,
      // and non-static references cannot be assignment targets.  Neither
      // can consts, so reference can't be replaced with const pointer.
      binaryRxnGenPair& rPair;
    
      void
      makeBinaryReactions(const leftContextClass& rNewContext,
			  const rightContextClass& rOtherContext) const
      {
	rPair.makeBinaryReactions(rNewContext,
				  rOtherContext);
      }

    public:

      leftRxnGenClass(rightFeatureClass& rRightFeature,
		      binaryRxnGenPair& rParentPair) :
	binaryRxnGen<leftContextClass, rightFeatureClass>(rRightFeature),
	rPair(rParentPair)
      {}
    };
    friend class leftRxnGenClass;
  
    leftRxnGenClass leftRxnGen;

    class rightRxnGenClass :
      public binaryRxnGen<rightContextClass, leftFeatureClass>
    {
      binaryRxnGenPair& rPair;

      void
      makeBinaryReactions(const rightContextClass& rNewContext,
			  const leftContextClass& rOtherContext) const
      {
	rPair.makeBinaryReactions(rOtherContext,
				  rNewContext);
      }

    public:

      rightRxnGenClass(leftFeatureClass& rLeftFeature,
		       binaryRxnGenPair& rParentPair) :
	binaryRxnGen<rightContextClass, leftFeatureClass>(rLeftFeature),
	rPair(rParentPair)
      {}
    };
    friend class rightRxnGenClass;
  
    rightRxnGenClass rightRxnGen;

  public:

    binaryRxnGenPair(leftFeatureClass& rLeftFeature,
		     rightFeatureClass& rRightFeature) :
      leftRxnGen(rRightFeature,
		 *this),
      rightRxnGen(rLeftFeature,
		  *this)
    {}

    virtual ~binaryRxnGenPair(void)
    {}

    // This virtual function must be supplied to tell how to construct
    // a reaction from the new contexts.
    virtual void
    makeBinaryReactions(const leftContextClass& rLeftContext,
			const rightContextClass& rRightContext) const = 0;

    /*! \brief Sub-generator accessors.

    These two functions allow adding the rxnGens to the features. */
    //@{
    typename leftFeatureClass::rxnGen*
    getLeftRxnGen()
    {
      return &leftRxnGen;
    }

    typename rightFeatureClass::rxnGen*
    getRightRxnGen()
    {
      return &rightRxnGen;
    }
    //@}
  };
}

#endif // BINARYRXNGEN_H

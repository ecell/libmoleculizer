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

#ifndef DIMERIZEEXTRAP_H
#define DIMERIZEEXTRAP_H

#include "dimer/dimerXcpt.hh"
#include "mzr/pchem.hh"

namespace dimer
{
  // Base class of rate extrapolators for dimerization reactions.
  class dimerizeExtrapolator
  {
  public:
    virtual
    ~dimerizeExtrapolator(void)
    {}
    
    virtual void
    setRate(bnd::siteParam leftParam,
	    bnd::siteParam rightParam,
	    double rate) = 0;
    
    virtual double
    getRate(const plx::siteInContext& rLeftContext,
	    const plx::siteInContext& rRightContext) const = 0;
  };

  // Dimerization rate extrapolator that does no extrapolation; uses the same
  // rate of dimerization for a given pair of binding site shapes, regardless
  // of dimerizing species.  Not recommended; included as an example only.
  class dimerizeNoExtrap :
    public dimerizeExtrapolator
  {
    typedef
    std::map<std::pair<bnd::siteParam, bnd::siteParam>, double>
    rateMapType;
    
    rateMapType rateMap;

  public:
    // Both for inserting default rates and for writing allosteric
    // rates over default rates.
    void
    setRate(bnd::siteParam leftParam,
	    bnd::siteParam rightParam,
	    double rate)
    {
      // For the time being, I'm storing the pair in only one order,
      // then searching for both orders when the pair is looked up.
      std::pair<rateMapType::iterator, bool> insertResult
	= rateMap.insert(std::make_pair(std::make_pair(leftParam,
						       rightParam),
					rate));
      if(! insertResult.second)
	{
	  insertResult.first->second = rate;
	}
    }

    double
    getRate(const plx::siteInContext& rLeftContext,
	    const plx::siteInContext& rRightContext) const
    {
      plx::cxSite cxLeft(rLeftContext);
      plx::cxSite cxRight(rRightContext);

      // Now we have to check for both orders in the pair of site shape
      // pointers.
      rateMapType::const_iterator iEntry
	= rateMap.find(std::make_pair(cxLeft.getSiteParam(),
				      cxRight.getSiteParam()));
      if(iEntry == rateMap.end())
	{
	  iEntry = rateMap.find(std::make_pair(cxRight.getSiteParam(),
					       cxLeft.getSiteParam()));
	  if(iEntry == rateMap.end())
	    throw missingDimerizeRateXcpt(cxLeft,
					  cxRight);
	}

      return iEntry->second;
    }
  };

  // Dimerization rate extrapolator that uses the masses of the
  // dimerizing species.
  class dimerizeMassExtrap :
    public dimerizeExtrapolator
  {
    typedef
    std::map<std::pair<bnd::siteParam, bnd::siteParam>, double>
    invMapType;

    double leftMass;
    double rightMass;

    invMapType invariantMap;
    
  public:
    // The masses given to this constructor are used to convert
    // to/from binding invariant to rate and back.
    dimerizeMassExtrap(double leftMolMass,
		       double rightMolMass) :
      leftMass(leftMolMass),
      rightMass(rightMolMass)
    {}

    // Both for inserting default rates and for writing allosteric
    // rates over default rates.
    // 
    // Storing the rate (binding invariant) with the key pair in only
    // one order; this means that search has to look for both orders.
    void
    setRate(bnd::siteParam leftParam,
	    bnd::siteParam rightParam,
	    double rate)
    {
      double invariant = mzr::bindingInvariant(rate,
					       leftMass,
					       rightMass);
      
      std::pair<invMapType::iterator, bool> insertResult
	= invariantMap.insert(std::make_pair(std::make_pair(leftParam,
						       rightParam),
					     invariant));
      if(! insertResult.second)
	{
	  insertResult.first->second = invariant;
	}
    }

    // For retrieving weight-corrected dimerization rates.
    double
    getRate(const plx::siteInContext& rLeftContext,
	    const plx::siteInContext& rRightContext) const
    {
      plx::cxSite cxLeft(rLeftContext);
      plx::cxSite cxRight(rRightContext);

      // Since the rates were stored with the key pair in only one
      // order, we have to check for both orders when looking it up.
      invMapType::const_iterator iEntry
	= invariantMap.find(std::make_pair(cxLeft.getSiteParam(),
					   cxRight.getSiteParam()));
      if(iEntry == invariantMap.end())
	{
	  iEntry = invariantMap.find(std::make_pair(cxRight.getSiteParam(),
						    cxLeft.getSiteParam()));
	  if(iEntry == invariantMap.end())
	    throw missingDimerizeInvariantXcpt(cxLeft,
					       cxRight);
	}

      return mzr::bindingRate(iEntry->second,
			      cxLeft.getPlexWeight(),
			      cxRight.getPlexWeight());
    }
  };
}

#endif // DIMERIZEEXTRAP_H

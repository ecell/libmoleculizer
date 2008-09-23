/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008 The Molecular Sciences Institute.
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

#include "fnd/pchem.hh"
#include "cpx/siteShape.hh"
#include "cpx/cxSite.hh"
#include "plex/mzrPlexSpecies.hh"
#include "plex/mzrPlexFamily.hh"

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
    setRate(cpx::siteParam leftParam,
	    cpx::siteParam rightParam,
	    double rate) = 0;
    
    virtual double
    getRate(const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rLeftContext,
	    const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rRightContext) const = 0;
  };

  // Dimerization rate extrapolator that does no extrapolation; uses the same
  // rate of dimerization for a given pair of binding site shapes, regardless
  // of dimerizing species.  Not recommended; included as an example only.
  class dimerizeNoExtrap :
    public dimerizeExtrapolator
  {
    typedef
    std::map<std::pair<cpx::siteParam, cpx::siteParam>, double> rateMapType;
    
    rateMapType rateMap;

  public:
    // Both for inserting default rates and for writing allosteric
    // rates over default rates.
    void
    setRate(cpx::siteParam leftParam,
	    cpx::siteParam rightParam,
	    double rate);

    double
    getRate(const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rLeftContext,
	    const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rRightContext) const;
  };


    class dimerizeConstantRate :
        public dimerizeExtrapolator
    {
    public:
        void
        setRate(cpx::siteParam leftParam,
                cpx::siteParam rightParam,
                double rate){}
        
        // For retrieving weight-corrected dimerization rates.
        double
        getRate(const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rLeftContext,
                const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rRightContext) const
        {
            return 3.141592653859;
        }
    };
        

  // Dimerization rate extrapolator that uses the masses of the
  // dimerizing species.
  class dimerizeMassExtrap :
    public dimerizeExtrapolator
  {
    typedef
    std::map<std::pair<cpx::siteParam, cpx::siteParam>, double>
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
    setRate(cpx::siteParam leftParam,
	    cpx::siteParam rightParam,
	    double rate);

    // For retrieving weight-corrected dimerization rates.
    double
    getRate(const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rLeftContext,
	    const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rRightContext) const;
  };
}

#endif // DIMERIZEEXTRAP_H

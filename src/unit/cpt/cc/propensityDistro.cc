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

#include <functional>
#include <algorithm>
#include "cpt/cptReaction.hh"
#include "cpt/propensityDistro.hh"

namespace cpt
{

  class badTargetPropensityXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(double targetPropensity,
	  double achievedPropensity)
    {
      std::ostringstream msgStream;
      msgStream << "Only got to propensity "
		<< achievedPropensity
		<< " but needed propensity "
		<< targetPropensity
		<< " getting sample of propensity distribution.";
      return msgStream.str();
    }
  public:
    badTargetPropensityXcpt(double targetPropensity,
			    double achievedPropensity) :
      utl::xcpt(mkMsg(targetPropensity,
		      achievedPropensity))
    {}
  };

  class achievesTgt :
    public std::unary_function<propensityDistro::value_type, bool>
  {
    const long double tgtProp;
    long double& rCumProp;
  public:
    achievesTgt(long double targetPropensity,
		long double& rCumulativePropensity) :
      tgtProp(targetPropensity),
      rCumProp(rCumulativePropensity)
    {}

    bool
    operator()(const argument_type& rPropRxnPair) const
    {
      rCumProp += rPropRxnPair.first;
      return (tgtProp <= rCumProp);
    }
  };
  
  propensityDistro::iterator
  propensityDistro::
  sample(void)
    throw(utl::xcpt)
  {
    long double tgtPropensity = uniform.sample() * getTotalPropensity();

    long double cumPropensity = 0.0;

    iterator iSample =
      std::find_if(begin(),
		   end(),
		   achievesTgt(tgtPropensity,
			       cumPropensity));

    if(end() == iSample) throw badTargetPropensityXcpt(tgtPropensity,
						       cumPropensity);
    return iSample;
  }

  void
  propensityDistro::
  addReaction(cptReaction* pCptReaction)
  {
    // Propensity 0.0 goes along with the initial totalPropensity of 0.0.
    push_front(value_type(0.0,
			  pCptReaction));
    pCptReaction->setDistroEntry(begin());
  }
}

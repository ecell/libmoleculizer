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

#ifndef DIMERXCPT_H
#define DIMERXCPT_H

#include "mzr/mzrXcpt.hh"
#include "plex/cxSiteParam.hh"
#include "plex/cxBindingParam.hh"

namespace dimer
{
  // Thrown by dimerizeNoExtrap extrapolator when the extrapolator is asked to
  // construct rates for site shape pairs that have not been given a nominal
  // rate.
  class missingDimerizeRateXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const plx::cxSite& cxLeft,
	  const plx::cxSite& cxRight);

  public:
    missingDimerizeRateXcpt(const plx::cxSite& cxLeft,
			    const plx::cxSite& cxRight) :
      mzr::mzrXcpt(mkMsg(cxLeft,
			 cxRight))
    {}
  };

  // Thrown by dimerizeMassExtrap extrapolator when the extrapolator is asked
  // to construct rates for site shape pairs that have not been given a
  // nominal rate.
  class missingDimerizeInvariantXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const plx::cxSite& cxLeft,
	  const plx::cxSite& cxRight);

  public:
    missingDimerizeInvariantXcpt(const plx::cxSite& cxLeft,
			 const plx::cxSite& cxRight) :
      mzr::mzrXcpt(mkMsg(cxLeft,
			 cxRight))
    {}
  };

  // Thrown by decomposeNoExtrap extrapolator when the extrapolator is asked
  // to construct a decomposition rate for site shape pairs that have not been
  // given a nominal rate.
  class missingDecomposeRateXcpt : public mzr::mzrXcpt
  {
    static std::string
    mkMsg(const plx::cxBinding& cx);

  public:
    missingDecomposeRateXcpt(const plx::cxBinding& cx) :
      mzr::mzrXcpt(mkMsg(cx))
    {}
  };
}

#endif // DIMERXCPT_H

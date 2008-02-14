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

#ifndef ODIEAPP_H
#define ODIEAPP_H

#include "utl/domJob.hh"
#include "rk4util/polymap.hh"
#include "odi/dump.hh"
#include "odi/odieGlue.hh"
#include "odi/odieParse.hh"

namespace odie
{
  class odieApp :
    public utl::dom::domBatchApp
  {
    // The state.  Alas, cannot use vector practically for concentrations.
    double* pConcentrations;
    double now;

    double dumpTimeDelta;
    double stopTime;

    // Control parameters for adaptive step size.
    double epsilonAbs;
    double epsilonRel;
    double stateScalar;
    double derivScalar;

    // Dumping.
    std::list<dumpStream> dumpStreams;
    void dumpHeaders(void) const;
    void dump(void) const;
    
    rk4util::polymap<double> derivative;

    // Construction of derivatives from parsed reactions.
    void makeDerivatives(std::vector<parserReaction>& rReactions);

  public:
    odieApp(int argc,
	    char** argv,
	    xmlpp::Document* pDoc) throw(std::exception);

    ~odieApp(void)
    {
      delete [] pConcentrations;
    }
    
    int
    run(void) throw(std::exception);
  };
}

#endif // ODIEAPP_H

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

#ifndef TAUAPP_H
#define TAUAPP_H

#include <sstream>
#include <vector>
#include <list>
#include "utl/domJob.hh"
#include "rk4util/rk4util.hh"
#include "rk4/rk4tau.hh"
#include "rk4/dump.hh"
#include "rk4/tauParse.hh"

namespace rk4tau
{
  class rejectionLimitXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(int rejectionLimit,
	  double when)
    {
      std::ostringstream msgStream;
      msgStream << "Rejected more than "
		<< rejectionLimit
		<< " Poisson samples at time "
		<< when
		<< ".";
      return msgStream.str();
    }
  public:
    rejectionLimitXcpt(int rejectionLimit,
		       double when) :
      utl::xcpt(mkMsg(rejectionLimit,
		      when))
    {}
  };

  class tauApp :
    public utl::dom::domBatchApp
  {
    // The maximum number of times Poisson samples can be rejected
    // due to getting a negative population.  A fatal error occurs
    // when this limit is reached.
    static const int rejectionLimit = 1000;
  
    double now;
    double dumpTime;
    double dumpInterval;
    double stopTime;

    // This is bound to the "absolute error" algorithm, which seems to be
    // the only reasonable choice among those I've implemented, and perhaps
    // the most reasonable choice at all.
    double epsilon;

    // This is for seeding the random number generator.  I will need my
    // handy-dandy hash function to convert this to an integer in an
    // appropriate way.
    std::string seedString;

    // Avogadro's number times the volume.
    //
    // I should be able to do away with this, along with all molar
    // conversions at simulation time, by just setting up the equations
    // in terms of numbers of molecules.
    double molarConst;

    // Right now, this is serving not only to receive the initial populations
    // from the parser, but also as a back-alley of communication between the
    // solving and dumping: the populations here are regularly updated, and the
    // dumpStreams access them here.
    std::vector<int> populations;

    polymap<double> derivatives;
    polymap<int> popsForCounts;

    // This has trouble being a vector, since vector push_back requires
    // an assignment operator, and dumpStream can't be assigned, due to
    // its containing refereneces.
    //
    // An alternative change would be to convert the references to
    // pointers.
    std::list<dumpStream> dumpStreams;

    void
    seedRandom(void) const;
  
    void
    dumpHeaders(void) const;

    void
    dump(void) const;

    // The step adaption algorithm. This is installed directly by the
    // parser, as is much of the private member data in this application
    // class.
    cashKarpAdaptor* pAdaptor;

    // Converts the reaction counts into molar reaction amounts,
    // then computes the Cash-Karp 4th order estimate of the gradient
    // over [t, t+delta], and the Cash-Karp error estimate, which
    // is put through the step-adaptor to get a new
    // value of delta, which is the returned value.
    double
    estimate(const std::vector<int>& rReactionCounts,
	     double delta,
	     std::vector<double>& rGradient);

    // Does a stochastic leap using the given estimated gradient.
    void
    leap(const std::vector<double>& rGradient,
	 double delta,
	 std::vector<int>& rReactionCounts, // in/out.
	 std::vector<int>& rPopulations, // out
	 int& rTotalReactionCount) throw(utl::xcpt); // out
  public:
    tauApp(int argc,
	   char** argv,
	   xmlpp::Document* pDoc) throw(std::exception);
    
    tauApp(void) :
      now(0.0),
      dumpTime(0.0),
      dumpInterval(0.0),
      stopTime(0.0),
      epsilon(1.0),
      seedString("rk4tau"),
      molarConst(1.0),
      pAdaptor(0)
    {}

    ~tauApp(void)
    {
      delete pAdaptor;
    }

    // Parsing routines.  Could I factor this out???
    void
    makeDerivatives(std::vector<parserReaction>& rReactions);

    void
    makeRateConverter(polymap<double>& rConverterPoly,
			      std::vector<parserReaction>& rReactions);
    void
    makeCountConverter(std::vector<parserReaction>& rReactions);

    void
    parseDomInput(xmlpp::Document* pDoc) throw(std::exception);

    int
    run(void) throw(std::exception);
  };
}

#endif // TAUAPP_H

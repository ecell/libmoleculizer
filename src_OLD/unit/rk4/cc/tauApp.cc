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
#include "sampleDist/sampleDist.hh"
#include "utl/message.hh"
#include "rk4util/linAlg.hh"
#include "rk4/tauApp.hh"

namespace rk4tau
{
  // This is duplicated from libMzr, which we can't link into this app.
  // An alternative would be to make
  // some sort of overall utilities library, possibly including domUtils?
  class linearHash
  {
    // These will need to be adjusted, I expect.  Or maybe not.
    static const size_t multiplier;
    static const size_t summand;
  public:
    size_t operator()(const size_t& rData) const;
    size_t operator()(const std::string& rString) const;
  };

  size_t
  linearHash::operator()(const size_t& rData) const
  {
    return (rData * multiplier) + summand;
  }
  

  // These will need to be adjusted, I expect.  Or maybe not.
  const size_t linearHash::multiplier = 2897564231ul;
  const size_t linearHash::summand = 3248630751ul;

  // This could be somewhat templatized....
  // It also appears to be defined in tauApp.cc.
  class charHashAccum : public std::unary_function<char, void>
  {
    size_t& rValue;
    linearHash lh;
  public:
    charHashAccum(size_t& rHashValue) :
      rValue(rHashValue)
    {
    }

    void
    operator()(char c) const
    {
      rValue = lh(rValue + lh((size_t) c));
    }
  };

  size_t
  linearHash::operator()(const std::string& rString) const
  {
    size_t hashValue = 42;
    std::for_each(rString.begin(),
		  rString.end(),
		  charHashAccum(hashValue));
    return hashValue;
  }

  void
  tauApp::seedRandom(void) const
  {
    // Seed the random number generator with the integer.
    linearHash lh;
    sampleDist::seedUniformSampler(lh(seedString));
  }

  class doStreamDump :
    public std::unary_function<dumpStream, void>
  {
    const std::vector<int>& rPops;
    double curTime;
  public:
    doStreamDump(const std::vector<int>& rSpeciesPops,
		 double dumpTime) :
      rPops(rSpeciesPops),
      curTime(dumpTime)
    {}

    void
    operator()(const dumpStream& rDumpStream) const
    {
      rDumpStream.dump(rPops,
		       curTime);
    }
  };

  void
  tauApp::dump(void) const
  {
    for_each(dumpStreams.begin(),
	     dumpStreams.end(),
	     doStreamDump(populations,
			  now));
  }

  void
  tauApp::dumpHeaders(void) const
  {
    for_each(dumpStreams.begin(),
	     dumpStreams.end(),
	     mem_fun_ref(&dumpStream::dumpHeaders));
  }

  // This could be scalar multiplication, if scalar multiplication were
  // flexible enough.  Alternatively, counts could be manipulated as
  // doubles.
  class countToConc :
    public std::unary_function<int, double>
  {
    double molarConst;
  public:
    countToConc(double molarConstant) :
      molarConst(molarConstant)
    {}

    double
    operator()(int count) const
    {
      return ((double) count) / molarConst;
    }
  };

  double
  tauApp::estimate(const std::vector<int>& rReactionCounts,
		   double delta,
		   std::vector<double>& rGradient)
  {
    int numReactions = rReactionCounts.size();
  
    // Compute reaction amounts from reaction counts; i.e. convert the
    // counts to concentrations.
    std::vector<double> reactionAmts(numReactions);
    transform(rReactionCounts.begin(),
	      rReactionCounts.end(),
	      reactionAmts.begin(),
	      countToConc(molarConst));

    // Invoke the core routine to compute a consensus derivative over
    // the next time inverval of length delta, along with the error
    // estimator Delta_1.
    std::vector<double> Delta_1(numReactions);
    cashKarpCore(reactionAmts,
		 derivatives,
		 delta,
		 rGradient,
		 Delta_1);

    // Use the "absolute error" step size adaptor.
    return pAdaptor->adaptDelta(Delta_1,
				reactionAmts,
				rGradient,
				delta);
  }

  class trackSamplePoisson :
    public std::unary_function<double, int>
  {
    int& rReactionCount;
  public:
    trackSamplePoisson(int& rTotalReactionCount) :
      rReactionCount(rTotalReactionCount)
    {}

    int
    operator()(double rate) const
    {
      int count = sampleDist::samplePoisson(rate);
      rReactionCount += count;
      return count;
    }
  };

  void
  tauApp::leap(const std::vector<double>& rGradient,
	       double delta,
	       std::vector<int>& rReactionCounts,
	       std::vector<int>& rNewPops,
	       int& rTotalReactionCount) throw(utl::xcpt)
  {
    int numReactions = rReactionCounts.size();

    // Convert the gradient (rates of increase of the molar reaction
    // amounts) to Poisson rates (rates of increase of the numbers
    // of reactions) by multiplying by the molar constant.
    std::vector<double>poissonRates(numReactions);
    scalarMult(rGradient,
	       delta * molarConst,
	       poissonRates);

    // Sample Poisson to get counts for each reaction over the next
    // delta time interval.  Reject the sample if executing that many
    // reactions would make any of the species' populations negative.
    bool popsAccepted = false;
    int rejectionCount = 0;
    while(! popsAccepted)
      {
	int newTotalReactionCount = rTotalReactionCount;

	// Sample Poisson to get the number of reactions of each type.
	std::vector<int> reactionCountDeltas(numReactions);
	transform(poissonRates.begin(),
		  poissonRates.end(),
		  reactionCountDeltas.begin(),
		  trackSamplePoisson(newTotalReactionCount));

	// Increment the reaction counts.
	std::vector<int> newReactionCounts(numReactions);
	transform(reactionCountDeltas.begin(),
		  reactionCountDeltas.end(),
		  rReactionCounts.begin(),
		  newReactionCounts.begin(),
		  std::plus<int>());

	// Convert the new reaction counts into populations.
	popsForCounts.evaluate(newReactionCounts,
			       rNewPops);

	// Test to see that all the populations are non-negative.
	// The Poisson sample is rejected if this is not the case.
	std::vector<int>::iterator iBadPop
	  = find_if(rNewPops.begin(),
		    rNewPops.end(),
		    std::bind2nd(std::less<int>(), 0));

	popsAccepted = (iBadPop == rNewPops.end());
	if(! popsAccepted)
	  {
	    if(rejectionLimit < rejectionCount++)
	      throw rejectionLimitXcpt(rejectionLimit,
				       now);
	  }
	else
	  {
	    rReactionCounts = newReactionCounts;
	    rTotalReactionCount = newTotalReactionCount;
	  }
      }
  }

  class runStatsMsg :
    public utl::message
  {
    static std::string
    mkMsg(int leaps,
	  int dumps,
	  int reactions)
    {
      std::ostringstream msgStream;
      msgStream << "Executed "
		<< leaps
		<< " leaps, "
		<< dumps
		<< " dumps, and "
		<< reactions
		<< " reactions.";
      return msgStream.str();
    }
  public:
    runStatsMsg(int leapCount,
		int dumpCount,
		int totalReactionCount) :
      utl::message(mkMsg(leapCount,
			    dumpCount,
			    totalReactionCount))
    {}
  };

  int
  tauApp::run(void) throw(std::exception)
  {
    int leapCount = 0;
    int dumpCount = 0;
    int totalReactionCount = 0;
    
    seedRandom();
    dumpHeaders();
    dump();
    dumpCount++;

    // This is unfortunately confusing.  This is the number of kinds of
    // reactions.
    int numReactions = derivatives.size();

    double delta = dumpInterval / 2.0;
    double nextDumpTime = dumpInterval;
    std::vector<int> reactionCounts(numReactions, 0);
    bool stillSimulating = true;
    while(stillSimulating)
      {
	bool stopNextDump = (stopTime <= nextDumpTime);
	if(stopNextDump) nextDumpTime = stopTime;
      
	std::vector<double> gradient(numReactions);
	// Reduce the size of delta until it looks good.
	// We may want to move forward when the estimate
	// for the next delta is only some fraction (<1)
	// of the current delta, instead of strictly larger
	// than the current delta.
	double newDelta = estimate(reactionCounts,
				   delta,
				   gradient);
	while(newDelta < delta)
	  {
// 	    std::cerr << "regressed; now = "
// 		      << now
// 		      << "; delta = "
// 		      << delta
// 		      << "; newDelta = "
// 		      << newDelta
// 		      << "."
// 		      << std::endl;
	    
	    delta = newDelta;
	    newDelta = estimate(reactionCounts,
				delta,
				gradient);
	  }

// 	std::cerr << "PROGRESSED; now = "
// 		  << now
// 		  << "."
// 		  << std::endl;

	// Would a leap of size delta take us past the end of the
	// simulation?
	if(nextDumpTime <= now + delta)
	  {
	    // Adjust the size of delta to hit the dump time at the end
	    // of the leap.
	    // 
	    // This can be a very short time interval.  If it is, then
	    // it makes a very bad delta: the error estimator loses
	    // significance and reports error 0, leading to all sorts of
	    // trouble.
	    //
	    // I'm trying leaving delta alone and starting over with it
	    // at the dump time.
	    double dumpDelta = nextDumpTime - now;

	    // Compute the consensus gradient for this (shorter) delta
	    // interval.
	    //
	    // We don't use the "next delta" from this estimation; since
	    // the interval dumpDelta can be quite small, the error
	    // estimator frequently loses significance.
	    estimate(reactionCounts,
		     dumpDelta,
		     gradient);

	    // Leap to the dump time.
	    leap(gradient,
		 dumpDelta,
		 reactionCounts,
		 populations,
		 totalReactionCount);
	    leapCount++;
	    now += dumpDelta;

	    // Do the dump.
	    dump();
	    dumpCount++;
	    nextDumpTime += dumpInterval;

	    // Restore the delta that we would have used if we had not
	    // stopped short for the dump.  This delta will be "vetted"
	    // as usual at the start of the next cycle.
	    //
	    // I'm hoping that it will work better than the delta
	    // derived from the short step, which suffers from loss of
	    // significance in the error estimator.
	    delta = newDelta;

	    stillSimulating = (! stopNextDump);
	  }
	else
	  {
	    leap(gradient,
		 delta,
		 reactionCounts,
		 populations,
		 totalReactionCount);
	    leapCount++;
	    now += delta;

	    // Adapt delta for the next cycle.
	    delta = newDelta;
	  }

	// We use the delta from the first estimation, in spite of the
	// fact that we may have taken a shorter leap to a dump point.
	// The error estimator tends to lose significance and
	// erroneously report "no error at all" when the leap interval
	// is extrememly short.
	delta = newDelta;
      }
    // Close all the dumpStreams.
    std::for_each(dumpStreams.begin(),
		  dumpStreams.end(),
		  mem_fun_ref(&dumpStream::close));

    runStatsMsg statsMsg(leapCount,
			 dumpCount,
			 totalReactionCount);
    statsMsg.issue();

    return 0;
  }
}

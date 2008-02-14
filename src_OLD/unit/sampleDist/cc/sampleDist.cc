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

#include <climits>
#include "sampleDist/sampleDist.hh"

namespace sampleDist
{
#if defined RANDOM_IS_RAND48

  void
  seedUniformSampler(int seed)
  {
    srand48(seed);
  }

#elif defined RANDOM_IS_RANDOM

  void
  seedUniformSampler(int seed)
  {
    srand((unsigned int) seed);
  }

#else
#error No uniform random number generator to sample.
#endif

  // Poisson sampler based on inverting the cumulative distribution
  // function.
  //
  // This works for small (<100) rates, but the exponential can
  // underflow for larger rates.
  int
  samplePoisson_1(double rate)
  {
    double uniform = sampleUniform();
    int count = 0;
    // There is a serious underflow problem here; for some rates that we
    // will encounter, say 1000, countProb becomes double-precision 0.
    double countProb = exp(- rate);
    double totProb = countProb;
    while(totProb < uniform)
      {
	countProb = (rate * countProb) / ((double) (++count));
	totProb += countProb;
      }
    return count;
  }

  // This works ok for a wide range of rates, but it is O(rate) and
  // can be fairly slow.
  int
  samplePoisson_2(double rate)
  {
    int count = 0;
    double time = sampleExp();
    while(time < rate)
      {
	count++;
	time += sampleExp();
      }
    return count;
  }

  // This works fine for reasonable (< 100) rates.  It is faster than
  // summing exponential variates, and maybe even faster than summing
  // the cumulative distribution function.  But it has the same
  // underflow problem as the cumulative distribution function
  // technique.
  int
  samplePoisson_3(double rate)
  {
    int count = 0;
    double limit = exp(-rate);
    double uniProd = sampleUniform();
    while(limit < uniProd)
      {
	count++;
	uniProd = uniProd * sampleUniform();
      }
    return count;
  }

  // gammaln from "Numerical Recipes" converted to double and prettied up.
  static const double
  gamCoef[6] =
    {
      76.18009172947146,
      -86.50532032941677,
      24.01409824083091,
      -1.231739572450155,
      0.1208650973866179e-2,
      -0.5395239384953e-5
    };

  double
  gammaln(double x)
  {
    double tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);

    double ser = 1.000000000190015;
    {
      double v = x;
      for(int ndx = 0;
	  ndx <= 5;
	  ndx++)
	{
	  ser += gamCoef[ndx] / ++v;
	}
    }

    return log(2.5066822746310005 * ser / x) - tmp;
  }

  // This is the rejection-based method out of "Numerical Recipes," with
  // its remaining mysterious spots.
  int
  samplePoisson(double rate)
  {
    // Use the simplest O(rate) method for small rates.
    if(rate < 20.0) return samplePoisson_3(rate);
    else
      {
	double sq = sqrt(2.0 * rate);
	double logRate = log(rate);
	double g = (rate * logRate) - gammaln(rate + 1.0);

	double em = 0.0;
	{
	  double t = 0.0;
	  do
	    {
	      // Sample the Cauchy distribution. (N.R. calls it "Lorentz"
	      // which I've never heard elsewhere.)  Reject negative
	      // points right away.
	      //
	      // Wouldn't it be better just to use uniform density on
	      // [0, pi/2] rather than this rejection thing?
	      double y = 0.0;
	      do
		{
		  // "Sample" the arctangent.
		  y = tan(3.14159 * sampleUniform());
		  // Adjusted to the shape (mean and variance) of
		  // Poisson.
		  em = (sq * y) + rate;
		}
	      while(em < 0.0);

	      // Now we are working with the histogram of the Poisson
	      // density.  Here, we determine the lhs of the histogram
	      // bar lying above the sample point em.
	      em = floor(em);

	      // Here, we compute the ratio of the arctangent density
	      // and the histogram density above em.
	      t = 0.9 * (1.0 + (y * y))
		* exp(em * logRate - gammaln(em + 1.0) - g);
	    }
	  while(t < sampleUniform());
	}
	return (int) em;
      }
  }
}

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

#ifndef SAMPLEDIST_H
#define SAMPLEDIST_H

/*! \file sample.hh
  \ingroup mzrGroup
  \brief Samplers for simple probability distributions. */

#include <cstdlib>
#include <cmath>
#include <climits>
#include "utl/platform.hh"

// Samplers for well-known distributions: uniform, exponential, discrete
// uniform, and Poisson.
namespace sampleDist
{
#if defined RANDOM_IS_RAND48

  inline double
  sampleUniform(void)
  {
    // drand48 is supposed to sample [0.0, 1.0), and we want to be able
    // to take the log of the return value of this function.
    return 1.0 - drand48();
  }

#elif defined RANDOM_IS_RANDOM

  inline double
  sampleUniform(void)
  {
    static const double dLongMax = (double) LONG_MAX;
    return ((double) random()) / dLongMax;
  }

#else
#error No uniform random number generator is known.
#endif

  void
  seedUniformSampler(int seed);


  /*! \ingroup mzrGroup
    \brief Samples the standard exponential distribution. */
  inline double
  sampleExp(void)
  {
    return - log(sampleUniform());
  }

  /*! \ingroup mzrGroup
    \brief Samples the uniform distribution on the integers {0,...,popSize}. */
  inline int
  samplePop(int popSize)
  {
#if defined RANDOM_IS_RAND48
    return lrand48() % popSize;
#elif defined RANDOM_IS_RANDOM
    return random() % popSize;
#else
#error No integer random number generator is known.
#endif
  }

  // Sample time is now O(rate); could be O(sqrt(rate)) with very little
  // elaboration.
  int
  samplePoisson(double rate);
}

#endif // SAMPLEDIST_H

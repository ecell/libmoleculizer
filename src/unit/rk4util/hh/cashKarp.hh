/////////////////////////////////////////////////////////////////////////////
// rk4 - a 4th order adaptive Runge Kutta solver for polynomial ODE's.
// Copyright (C) 2002  Walter Lawrence (Larry) Lok.
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

#ifndef CASHKARP_H
#define CASHKARP_H

#include "rk4util/polymap.hh"

namespace rk4util
{
  // This is the core routine for both rk4tau and rk4ode.
  void
  cashKarpCore(const std::vector<double>& rValues,
	       const polymap<double>& rDerivs,
	       double delta,
	       std::vector<double>& rGradient,
	       std::vector<double>& rDelta_1);

  // Used in recomputing delta in each cycle to get the minimum
  // of the error ratios over all the variables.
  double
  minRatio(const std::vector<double>& rBigDelta_0,
	   const std::vector<double>& rBigDelta_1);

  // Base class for rk4 step adaptor, which computes the length
  // of the next interval from the error estimate given by the
  // core routine.
  class cashKarpAdaptor
  {
  public:
    virtual ~cashKarpAdaptor(void)
    {}
  
    virtual double
    adaptDelta(const std::vector<double>& rDelta_1,
	       const std::vector<double>& rValues,
	       const std::vector<double>& rDerivs,
	       double delta) const = 0;
  };

  // The "absolute error" adaptor.
  class ckStepAdaptor :
    public cashKarpAdaptor
  {
    // The maximum ratio (here 10) by which delta may expand, raised to
    // the power 1/"expansion power".
    static const double maxExpansion;
    static const double expansionPower;

    double BigDelta_0;
  
  public:
    ckStepAdaptor(double stepError) :
      BigDelta_0(stepError)
    {}

    virtual double
    adaptDelta(const std::vector<double>& rBigDelta_1,
	       const std::vector<double>& rValues,
	       const std::vector<double>& rDerivs,
	       double delta) const;
  };

  // The "relative error" adaptor. (Seems like this one was catastrophic.)
  class ckRelErrAdaptor :
    public cashKarpAdaptor
  {
    class makeDelta_0;

    double eps;

  public:
    ckRelErrAdaptor(const double& epsilon) :
      eps(epsilon)
    {}

    virtual double
    adaptDelta(const std::vector<double>& rBigDelta_1,
	       const std::vector<double>& rValues,
	       const std::vector<double>& rDerivs,
	       double delta) const;
  };

  // A form of relative error adaptor; recommended by "Numerical
  // Recipes" for stiff equations.  But since reaction amounts increase
  // constantly, and the error in the populations must not grow
  // continuously, it seems as though any form of relative error
  // estimation in the reaction amounts is rather useless.
  class ckStiffAdaptor :
    public cashKarpAdaptor
  {
    class makeDelta_0;

    // This is the constant mentioned on p.737,
    double BigC;
    double eps;
  public:
    ckStiffAdaptor(double Cthreshold,
		   double epsilon) :
      BigC(Cthreshold),
      eps(epsilon)
    {}

    virtual double
    adaptDelta(const std::vector<double>& rBigDelta_1,
	       const std::vector<double>& rValues,
	       const std::vector<double>& rDerivs,
	       double delta) const;
  };
}

#endif // CASHKARP_H

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

#include <cmath>
#include <float.h>
#include <algorithm>
#include <functional>
#include "fnd/physConst.hh"
#include "rk4util/rk4util.hh"

namespace rk4util
{
  // Used in minRatio.
  static const double maxExpansion = 1.0e5;
  static const double expansionPower = 0.2;

  // These ai aren't useful for us, since we are solving autonomous
  // systems only.
  static const double ck_a_2 = 1.0/5.0;
  static const double ck_a_3 = 3.0/10.0;
  static const double ck_a_4 = 3.0/5.0;
  static const double ck_a_5 = 1.0;
  static const double ck_a_6 = 7.0/8.0;

  static const double ck_b_2_1 = 1.0/5.0;
  static const double ck_b_3_1 = 3.0/40.0;
  static const double ck_b_3_2 = 9.0/40.0;
  static const double ck_b_4_1 = 3.0/10.0;
  static const double ck_b_4_2 = - 9.0/10.0;
  static const double ck_b_4_3 = 6.0/5.0;
  static const double ck_b_5_1 = - 11.0/54.0;
  static const double ck_b_5_2 = 5.0/2.0;
  static const double ck_b_5_3 = - 70.0/27.0;
  static const double ck_b_5_4 = 35.0/27.0;
  static const double ck_b_6_1 = 1631.0/55296.0;
  static const double ck_b_6_2 = 175.0/512.0;
  static const double ck_b_6_3 = 575.0/13824.0;
  static const double ck_b_6_4 = 44275.0/110592.0;
  static const double ck_b_6_5 = 253.0/4096.0;

  static const double ck_c_1 = 37.0/378.0;
  static const double ck_c_2 = 0.0;
  static const double ck_c_3 = 250.0/621.0;
  static const double ck_c_4 = 125.0/594.0;
  static const double ck_c_5 = 0.0;
  static const double ck_c_6 = 512.0/1771.0;

  static const double ck_cs_1 = 2825.0/27648.0;
  static const double ck_cs_2 = 0.0;
  static const double ck_cs_3 = 18575.0/48384.0;
  static const double ck_cs_4 = 13525.0/55296.0;
  static const double ck_cs_5 = 277.0/14336.0;
  static const double ck_cs_6 = 1.0/4.0;

  // I've decided not to emit the value of the derivative at the current
  // solution point, even though it is produced as a by-product.  I had
  // assumed that the step-adaption algorithm might use it, but now
  // would just use the consensus derivative there instead.  But my most
  // likely step adaption algorithm only uses BigDelta_1:
  //
  // BigDelta_1 = (delta)(rGrad - rAltGrad).
  void
  cashKarpCore(const std::vector<double>& rValues,
	       const polymap<double>& rDerivs,
	       double delta,
	       std::vector<double>& rGrad,
	       std::vector<double>& rDelta_1)
  {
    int dimension = rValues.size();

    // Compute the derivatives at the current solution point.
    std::vector<double> d1(dimension);
    rDerivs.evaluate(rValues,
		     d1);

    // Compute k1, which I will call an "Euler point," as the calculation
    // is like a step in Euler's method.
    std::vector<double> k1(dimension);
    scalarMult(d1,
	       delta,
	       k1);

    // Compute d2, another derivative.
    std::vector<double> d2(dimension);
    {
      std::vector<double> vars2(dimension);
      {
	std::vector<double> sk1(dimension);
	scalarMult(k1, ck_b_2_1, sk1);
	vectorAdd(rValues, sk1, vars2);
      }
      rDerivs.evaluate(vars2,
		       d2);
    }

    // Compute k2, the next "Euler point."
    std::vector<double> k2(dimension);
    scalarMult(d2,
	       delta,
	       k2);

    // Compute the next derivative, d3.
    std::vector<double> d3(dimension);
    {
      std::vector<double> vars3_2(dimension);
      {
	std::vector<double> vars3_1(dimension);
	{
	  std::vector<double> sk1(dimension);
	  scalarMult(k1, ck_b_3_1, sk1);
	  vectorAdd(rValues, sk1, vars3_1);
	}
	std::vector<double> sk2(dimension);
	scalarMult(k2, ck_b_3_2, sk2);
	vectorAdd(vars3_1, sk2, vars3_2);
      }
      rDerivs.evaluate(vars3_2,
		       d3);
    }

    // Compute the next "Euler point" k3.
    std::vector<double> k3(dimension);
    scalarMult(d3,
	       delta,
	       k3);

    // Compute the next derivative, d4.
    std::vector<double> d4(dimension);
    {
      std::vector<double> vars4_3(dimension);
      {
	std::vector<double> vars4_2(dimension);
	{
	  std::vector<double> vars4_1(dimension);
	  {
	    std::vector<double> sk1(dimension);
	    scalarMult(k1, ck_b_4_1, sk1);
	    vectorAdd(rValues, sk1, vars4_1);
	  }
	  std::vector<double> sk2(dimension);
	  scalarMult(k2, ck_b_4_2, sk2);
	  vectorAdd(vars4_1, sk2, vars4_2);
	}
	std::vector<double> sk3(dimension);
	scalarMult(k3, ck_b_4_3, sk3);
	vectorAdd(vars4_2, sk3, vars4_3);
      }
      rDerivs.evaluate(vars4_3,
		       d4);
    }
  

    // Compute the next "Euler point," k4.
    std::vector<double> k4(dimension);
    scalarMult(d4,
	       delta,
	       k4);

    // Compute the next derivative, d5.
    std::vector<double> d5(dimension);
    {
      std::vector<double> vars5_4(dimension);
      {
	std::vector<double> vars5_3(dimension);
	{
	  std::vector<double> vars5_2(dimension);
	  {
	    std::vector<double> vars5_1(dimension);
	    {
	      std::vector<double> sk1(dimension);
	      scalarMult(k1, ck_b_5_1, sk1);
	      vectorAdd(rValues, sk1, vars5_1);
	    }
	    std::vector<double> sk2(dimension);
	    scalarMult(k2, ck_b_5_2, sk2);
	    vectorAdd(vars5_1, sk2, vars5_2);
	  }
	  std::vector<double> sk3(dimension);
	  scalarMult(k3, ck_b_5_3, sk3);
	  vectorAdd(vars5_2, sk3, vars5_3);
	}
	std::vector<double> sk4(dimension);
	scalarMult(k4, ck_b_5_4, sk4);
	vectorAdd(vars5_3, sk4, vars5_4);
      }
      rDerivs.evaluate(vars5_4,
		       d5);
    }

    // Compute the next "Euler point," k5.
    std::vector<double> k5(dimension);
    scalarMult(d5,
	       delta,
	       k5);

    // Compute the last derivative, d6.
    std::vector<double> d6(dimension);
    {
      std::vector<double> vars6_5(dimension);
      {
	std::vector<double> vars6_4(dimension);
	{
	  std::vector<double> vars6_3(dimension);
	  {
	    std::vector<double> vars6_2(dimension);
	    {
	      std::vector<double> vars6_1(dimension);
	      {
		std::vector<double> sk1(dimension);
		scalarMult(k1, ck_b_6_1, sk1);
		vectorAdd(rValues, sk1, vars6_1);
	      }
	      std::vector<double> sk2(dimension);
	      scalarMult(k2, ck_b_6_2, sk2);
	      vectorAdd(vars6_1, sk2, vars6_2);
	    }
	    std::vector<double> sk3(dimension);
	    scalarMult(k3, ck_b_6_3, sk3);
	    vectorAdd(vars6_2, sk3, vars6_3);
	  }
	  std::vector<double> sk4(dimension);
	  scalarMult(k4, ck_b_6_4, sk4);
	  vectorAdd(vars6_3, sk4, vars6_4);
	}
	std::vector<double> sk5(dimension);
	scalarMult(k5, ck_b_6_5, sk5);
	vectorAdd(vars6_4, sk5, vars6_5);
      }
      rDerivs.evaluate(vars6_5,
		       d6);
    }

    // In this way of doing things, there's no reason to compute
    // k6, the final "Euler point."

    // Construct the consensus derivative.
    {
      std::vector<double> cd5(dimension);
      {
	std::vector<double> cd4(dimension);
	{
	  std::vector<double> cd3(dimension);
	  {
	    std::vector<double> cd2(dimension);
	    {
	      std::vector<double> cd1(dimension);
	      scalarMult(d1, ck_c_1, cd1);

	      std::vector<double> sd2(dimension);
	      scalarMult(d2, ck_c_2, sd2);
	      vectorAdd(cd1, sd2, cd2);
	    }
	    std::vector<double> sd3(dimension);
	    scalarMult(d3, ck_c_3, sd3);
	    vectorAdd(cd2, sd3, cd3);
	  }
	  std::vector<double> sd4(dimension);
	  scalarMult(d4, ck_c_4, sd4);
	  vectorAdd(cd3, sd4, cd4);
	}
	std::vector<double> sd5(dimension);
	scalarMult(d5, ck_c_5, sd5);
	vectorAdd(cd4, sd5, cd5);
      }
      std::vector<double> sd6(dimension);
      scalarMult(d6, ck_c_6, sd6);
      vectorAdd(cd5, sd6, rGrad);
    }

    // Construct the alternative consensus derivative.
    std::vector<double> altGrad(dimension);
    {
      std::vector<double> cd5(dimension);
      {
	std::vector<double> cd4(dimension);
	{
	  std::vector<double> cd3(dimension);
	  {
	    std::vector<double> cd2(dimension);
	    {
	      std::vector<double> cd1(dimension);
	      scalarMult(d1, ck_cs_1, cd1);

	      std::vector<double> sd2(dimension);
	      scalarMult(d2, ck_cs_2, sd2);
	      vectorAdd(cd1, sd2, cd2);
	    }
	    std::vector<double> sd3(dimension);
	    scalarMult(d3, ck_cs_3, sd3);
	    vectorAdd(cd2, sd3, cd3);
	  }
	  std::vector<double> sd4(dimension);
	  scalarMult(d4, ck_cs_4, sd4);
	  vectorAdd(cd3, sd4, cd4);
	}
	std::vector<double> sd5(dimension);
	scalarMult(d5, ck_cs_5, sd5);
	vectorAdd(cd4, sd5, cd5);
      }
      std::vector<double> sd6(dimension);
      scalarMult(d6, ck_cs_6, sd6);
      vectorAdd(cd5, sd6, altGrad);
    }

    // Construct Delta_1.  Here, I'm using the fact that the y_n
    // cancel in the formula, so that rDelta_1 is just the difference
    // of the two "consensus" derivatives times delta.
    std::vector<double> derivDifference(dimension);
    {
      std::vector<double> negAltGrad(dimension);
      scalarMult(altGrad, -1.0, negAltGrad);
      vectorAdd(rGrad, negAltGrad, derivDifference);
    }
    scalarMult(derivDifference, delta, rDelta_1);
  }

  double
  minRatio(const std::vector<double>& rBigDelta_0,
	   const std::vector<double>& rBigDelta_1)
  {
    // This puts a limit on the factor by which step size can increase.
    double minRatio = maxExpansion;

    int varbNdx = (int) rBigDelta_1.size();
    while(0 < varbNdx--)
      {
	double D_0 = rBigDelta_0[varbNdx];
	if(D_0 < 0.0) D_0 = - D_0;

	double D_1 = rBigDelta_1[varbNdx];
	if(D_1 < 0.0) D_1 = - D_1;

	// As it turns out bigDelta_1 can come out 0, which would
	// mean that the error was very small.  This inequality
	// only rules out cases whose ratio would be > MAX_EXPANSION anyway,
	// and so would not affect the minRatio.
	if(D_0 < (D_1 * maxExpansion))
	  {
	    double errorRatio = D_0 / D_1;
	    if(errorRatio < minRatio) minRatio = errorRatio;
	  }
      }
    return pow(minRatio, expansionPower);
  }

  double
  ckStepAdaptor::adaptDelta(const std::vector<double>& rBigDelta_1,
			    const std::vector<double>& rValues,
			    const std::vector<double>& rDerivs,
			    double delta) const
  {
    double minimumRatio = maxExpansion;

    int varbNdx = (int) rBigDelta_1.size();
    while(0 < varbNdx--)
      {
	double D_1 = rBigDelta_1[varbNdx];
	if(D_1 < 0.0) D_1 = - D_1;

	if(BigDelta_0 < (D_1 * maxExpansion))
	  {
	    double errorRatio = BigDelta_0 / D_1;
	    if(errorRatio < minimumRatio) minimumRatio = errorRatio;
	  }
      }

    return delta * pow(minimumRatio, expansionPower);
  }

  class ckRelErrAdaptor::makeDelta_0 :
    public std::binary_function<double, double, double>
  {
    double epsilon;
    double delta;
    
  public:
    makeDelta_0(const double& rEpsilon,
		const double& rDelta) :
      epsilon(rEpsilon),
      delta(rDelta)
    {}

    double
    operator()(double value,
	       double deriv) const
    {
      // Get absolute values. Assuming delta is positive.
      if(0.0 <= value) value = -value;
      if(0.0 <= deriv) deriv = - deriv;

      return epsilon * (value + (delta * deriv));
    }
  };

  double
  ckRelErrAdaptor::adaptDelta(const std::vector<double>& rBigDelta_1,
			      const std::vector<double>& rValues,
			      const std::vector<double>& rDerivs,
			      double delta) const
  {
    std::vector<double> BigDelta_0;

    // Compute BigDelta_0 from the values and derivatives.
    transform(rValues.begin(),
	      rValues.end(),
	      rDerivs.begin(),
	      back_inserter(BigDelta_0),
	      makeDelta_0(eps,
			  delta));

    // Compute the rato.
    double ratio = minRatio(BigDelta_0,
			    rBigDelta_1);

    return delta * ratio;
  }

  class ckStiffAdaptor::makeDelta_0 :
    public std::unary_function<double, double>
  {
    double BigC;
    double eps;
    
  public:
    makeDelta_0(double Cthreshold,
		double epsilon) :
      BigC(Cthreshold),
      eps(epsilon)
    {}

    double
    operator()(const double& rValue) const
    {
      double val = rValue;
    
      if(val < 0.0) val = -val;
      if(val < BigC)
	{
	  val = BigC;
	}
      return eps * val;
    }
  };

  double
  ckStiffAdaptor::adaptDelta(const std::vector<double>& rBigDelta_1,
			     const std::vector<double>& rValues,
			     const std::vector<double>& rDerivs,
			     double curDelta) const
  {
    std::vector<double> BigDelta_0(rBigDelta_1.size());

    // Compute BigDelta_0 from the values and derivatives.
    transform(rValues.begin(),
	      rValues.end(),
	      BigDelta_0.begin(),
	      makeDelta_0(BigC,
			  eps));

    // Update delta.
    return curDelta * minRatio(BigDelta_0,
			       rBigDelta_1);
  }

  const double ckStepAdaptor::maxExpansion = 1.0e5;
  const double ckStepAdaptor::expansionPower = 0.2;
}

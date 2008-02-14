#include <iostream>
#include <complex>
#include "rk4util/polymap.hh"
#include "rk4util/linAlg.hh"

std::complex<double>
operator*(const int& rLeft,
	  const std::complex<double>& rRight)
{
  return std::complex<double>(rLeft, 0.0) * rRight;
}

std::complex<double>
operator*(const std::complex<double>& rRight,
	  const int& rLeft)
{
  return std::complex<double>(rLeft, 0.0) * rRight;
}


int main(int argc, char* argv[])
{
  // Construct the polynomial 3x^3 + 2xy + 1.
  rk4util::monomial xCubed;
  xCubed.setExponent(0, 3);
  rk4util::monomial xy;
  xy.setExponent(0, 1);
  xy.setExponent(1, 1);
  rk4util::monomial one;
  
  rk4util::polynomial<double> testPoly;
  testPoly.setCoeff(xCubed, 3.0);
  testPoly.setCoeff(xy, 2.0);
  testPoly.setCoeff(one, 1.0);

  // Evaluate testPoly at (3, 7)
  std::vector<double> testArg(2);
  testArg[0] = 3.0;
  testArg[1] = 7.0;

  std::cerr << "Test polynomial evaluates to "
	    << testPoly.evaluate(testArg)
	    << "."
	    << std::endl;

  // Extend scalars of testPoly to complex.
  rk4util::polynomial<std::complex<double> > testPolyCx;
  testPoly.extendScalars(testPolyCx);
  std::vector<std::complex<double> > testArgCx(testArg.size());
  rk4util::extendScalars(testArg,
			 testArgCx);

  std::cerr << "Test polynomial (complexified) evaluates to "
	    << testPolyCx.evaluate(testArgCx)
	    << "."
	    << std::endl;

  // Multiply complexified testPoly by itself and evaluate.
  rk4util::polynomial<std::complex<double> > testPolySquared;
  testPolyCx.multiply(testPolyCx,
		      testPolySquared);

  std::cerr << "Test polynomial(complexified) squared evaluates to "
	    << testPolySquared.evaluate(testArgCx)
	    << "."
	    << std::endl;

  // Add squared test polynomial to itself and evaluate.
  rk4util::polynomial<std::complex<double> > twiceTestPolySquared
    = testPolySquared.add(testPolySquared);
  
  std::cerr << "Twice test polynomial...squared evaluates to "
	    << twiceTestPolySquared.evaluate(testArgCx)
	    << "."
	    << std::endl;

  // Form polymap of all three polynomials.
  rk4util::polymap<std::complex<double> > testPm;
  testPm.push_back(testPolyCx);
  testPm.push_back(testPolySquared);
  testPm.push_back(twiceTestPolySquared);

  std::vector<std::complex<double> > polyMapResult
    = testPm.evaluate(testArgCx);

  std::cerr << "Poly map application result: ["
	    << polyMapResult[0]
	    << ", "
	    << polyMapResult[1]
	    << ", "
	    << polyMapResult[2]
	    << "]."
	    << std::endl;

  // Perform partial derivative wrt x of test poly and evaluate at (3, 7).
  rk4util::polynomial<double> pdx
    = testPoly.partialDerivative<double>(0);

  std::cerr << "Partial derivative wrt x result: "
	    << pdx.evaluate(testArg)
	    << std::endl;

  // Perform partial derivative wrt y of test poly and evaluate at (3, 7).
  rk4util::polynomial<double> pdy
    = testPoly.partialDerivative<double>(1);

  std::cerr << "Partial derivative wrt y result: "
	    << pdy.evaluate(testArg)
	    << std::endl;

  // What happens when you take a partial derivative wrt an off-the-wall
  // variable?
  rk4util::polynomial<double> pdoff
    = testPoly.partialDerivative<double>(42);

  std::cerr << "Partial derivative wrt off-the-wall result: "
	    << pdoff.evaluate(testArg)
	    << std::endl;

  // Testing gradient of a polynomial.
  // Note that the gradient polymap has to have the right number
  // of components (the number of variables.)
  rk4util::polymap<double> grad(2);
  testPoly.gradient(2,
		    grad);

  std::vector<double> gradEval
    = grad.evaluate(testArg);

  std::cerr << "Gradient application result: ["
	    << gradEval[0]
	    << ", "
	    << gradEval[1]
	    << "]."
	    << std::endl;

  // Testing Jacobian of a polymap.
  std::vector<rk4util::polymap<std::complex<double> > > jcb
    = testPm.jacobian<std::complex<double> >(2);

  std::vector<std::complex<double> > jcbEval0
    = jcb[0].evaluate(testArgCx);

  std::cerr << "Jacobian application result: ["
	    << jcbEval0[0]
	    << ", "
	    << jcbEval0[1]
	    << "]"
	    << std::endl;

  std::vector<std::complex<double> > jcbEval1
    = jcb[1].evaluate(testArgCx);

  std::cerr << "Jacobian application result: ["
	    << jcbEval1[0]
	    << ", "
	    << jcbEval1[1]
	    << "]"
	    << std::endl;

  std::vector<std::complex<double> > jcbEval2
    = jcb[2].evaluate(testArgCx);

  std::cerr << "Jacobian application result: ["
	    << jcbEval2[0]
	    << ", "
	    << jcbEval2[1]
	    << "]"
	    << std::endl;

  rk4util::polymap<std::complex<double> > grad2(2);
  testPolySquared.gradient(2,
			   grad2);
  std::vector<std::complex<double> > grad2Eval
    = grad2.evaluate(testArgCx);

  std::cerr << "Checking second row: ["
	    << grad2Eval[0]
	    << ", "
	    << grad2Eval[1]
	    << "]"
	    << std::endl;

  return 0;
}

  

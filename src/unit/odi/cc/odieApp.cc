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

#include "odi/odieApp.hh"

namespace odie
{
  int
  odieApp::run(void) throw(std::exception)
  {
    // Compute the Jacobian matrix for use by the implicit solver.
    std::vector<rk4util::polymap<double> > jacobian
      = derivative.jacobian<double>(derivative.size());

    // Package functional information up for gsl.
    odieParams params(derivative,
		      jacobian);

    // Set up the equations for gsl.
    gsl_odeiv_system system;
    system.function = &evalDeriv;
    system.jacobian = &evalJacobian;
    system.dimension = derivative.size();
    system.params = (void *) &params;

    // Create a Bulirsch-Stoer stepper.
    gsl_odeiv_step* pStep = gsl_odeiv_step_alloc(gsl_odeiv_step_bsimp,
						 derivative.size());

    // Create an implicit 4th order rk stepper.
    //     gsl_odeiv_step* pStep = gsl_odeiv_step_alloc(gsl_odeiv_step_rk4imp,
    //  						 derivative.size());

    // Create one of the other implicit steppers.
    //     gsl_odeiv_step* pStep = gsl_odeiv_step_alloc(gsl_odeiv_step_gear2,
    //   						 derivative.size());

    // Create adaptive step size controller.
    gsl_odeiv_control* pControl
      = gsl_odeiv_control_standard_new(epsilonAbs,
				       epsilonRel,
				       stateScalar,
				       derivScalar);

    // Create "evolver" that runs the stepper and step size controller
    // to move forward in time.
    gsl_odeiv_evolve* pEvolver
      = gsl_odeiv_evolve_alloc(derivative.size());

    // Dump headers to the output files.
    dumpHeaders();

    double stepSize = 1.0;
    int status = GSL_SUCCESS;
    // I intend at this point for dumpTimeDelta to default to -1.0
    // unless overridden by an optional, positive, user-provided value.
    if(0.0 < dumpTimeDelta)
      {
	// Run the solver from dump time to dump time.
	double nextDumpTime = now + dumpTimeDelta;
	while((GSL_SUCCESS == status)
	      && (now < stopTime))
	  {
	    dump();
	    
	    // Don't run over the end of the simulation looking
	    // for the end of the next dump cycle!
	    double cycleEndTime = nextDumpTime;
	    if(stopTime < cycleEndTime)
	      {
		cycleEndTime = stopTime;
	      }

	    // Run the simulation until it reaches the end of the
	    // dump cycle, or the end of the simulation, whichever
	    // comes first.
	    while((GSL_SUCCESS == status)
		  && (now < cycleEndTime))
	      {
		status
		  = gsl_odeiv_evolve_apply(pEvolver,
					   pControl,
					   pStep,
					   &system,
					   &now,
					   cycleEndTime,
					   &stepSize,
					   pConcentrations);
	      }

	    nextDumpTime += dumpTimeDelta;
	  }
      }
    else
      {
	// Let the solver decide when to dump.
	while((GSL_SUCCESS == status)
	      && (now < stopTime))
	  {
	    dump();

	    status
	      = gsl_odeiv_evolve_apply(pEvolver,
				       pControl,
				       pStep,
				       &system,
				       &now,
				       stopTime,
				       &stepSize,
				       pConcentrations);
	  }
      }
    return 0;
  }

  class doStreamDump :
    public std::unary_function<dumpStream, void>
  {
    const double* pConc;
    double now;
  public:
    doStreamDump(const double* pConcentrations,
		 double curTime) :
      pConc(pConcentrations),
      now(curTime)
    {
    }

    void
    operator()(const dumpStream& rDumpStream) const
    {
      rDumpStream.dump(pConc,
		       now);
    }
  };

  void
  odieApp::dump(void) const
  {
    for_each(dumpStreams.begin(),
	     dumpStreams.end(),
	     doStreamDump(pConcentrations,
			  now));
  }

  void
  odieApp::dumpHeaders(void) const
  {
    for_each(dumpStreams.begin(),
	     dumpStreams.end(),
	     std::mem_fun_ref(&dumpStream::dumpHeaders));
  }
}

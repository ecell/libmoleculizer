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

#ifndef PARAMETRIZER_H
#define PARAMETRIZER_H

#include "utl/dom.hh"
#include "mzr/moleculizer.hh"

namespace mzr
{
  class parametrizer :
    public moleculizer
  {
  public:

    // The point of this program is to reparametrize a reaction network, in
    // the form of a moleculizer-state dump, using new rates etc.  substituted
    // into the moleculizer-input from which the original simulation was run.
    // 
    // It parses moleculizer-input as usual, then parses tagged species from
    // moleculizer-state.  The populations on the tagged species are
    // essentially ignored, and the populations of explicit species in the
    // moleculizer-input are treated as usual, so the result is the same
    // simulation as before, but with an expanded set of species and
    // reactions.
    parametrizer(int argc,
		 char** argv,
		 xmlpp::Document* pMoleculizerInput,
		 xmlpp::Document* pMoleculizerState)
      throw(std::exception);

    ~parametrizer(void);

    // This basically does nothing except "makeDomOutput" to produce the
    // reparametrized version of the input moleculizer-state file.  Does not
    // run the moleculizer event queue.
    int
    run(void) throw(std::exception);
  };
}

#endif

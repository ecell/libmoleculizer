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

#ifndef PHYSCONST_H
#define PHYSCONST_H

/*! \file physConst.hh
  \ingroup chemGroup
  \brief Compiled-in physical constants.

  I've generally tried to avoid these, preferring to put in such things
  in scripts.  

  Sometimes, a compiled-in placeholder is necessary.
  In the alpha units, I have a header alphaConst.hh that gives a number of
  such placeholders, but they are not intended for use.

  \todo At present there is no warning if a placeholder physical constant
  has not been filled in.  How to do it?  Put check in prepareToRun?

  Avogadro's number is needed for basic reaction rate calculations
  in many contexts, and it is not likely to need revision. */

#define AVOGADROS_NUMBER (6.02e23)

#endif

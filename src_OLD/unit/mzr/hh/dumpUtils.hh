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

#ifndef DUMPUTILS_H
#define DUMPUTILS_H

#include "utl/dom.hh"

namespace mzr
{
  // This may want to go into utl; I've never really "utilitized" output.

  // Inserts a stereotyped double-valued parameter element, with additional
  // elements giving the parameter value in scientific notation
  //
  // This is basically part of a fix to introduce scientific notation for all
  // parameter values so as to be able, using XSLT, to generate SBML and other
  // formats that have the fraction and exponent in separate XML constructs.
  void
  addDoubleParamChild(xmlpp::Node* pParentNode,
		      const std::string& rChildName,
		      const std::string& rParameterName,
		      double parameterValue);
}

#endif // DUMPUTILS_H

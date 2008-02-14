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

#ifndef CPT_SIMTIMEDUMPABLE_H
#define CPT_SIMTIMEDUMPABLE_H

#include "fnd/dumpable.hh"

namespace cpt
{
  class cptApp;
  
  /*! \ingroup dumpGroup
    \brief Dumps the simulation time.

    Automatically used as the first field by most dump commands. */
  class simTimeDumpable :
    public fnd::dumpable<fnd::basicDumpable::dumpArg>
  {
    cptApp& rApp;
    
  public:
    simTimeDumpable(cptApp& rCptApp);
      
    virtual void
    doDump(const fnd::basicDumpable::dumpArg& rDumpArg) const;
  };
}

#endif // CPT_SIMTIMEDUMPABLE_H

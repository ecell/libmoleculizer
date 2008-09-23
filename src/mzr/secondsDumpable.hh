/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008 The Molecular Sciences Institute.
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef MZR_SECONDSDUMPABLE_H
#define MZR_SECONDSDUMPABLE_H

#include "fnd/dumpable.hh"

namespace mzr
{
  class mzrUnit;
  
  /*! \ingroup dumpGroup
    \brief Dumps dumps clock time since program start in whole seconds.

    This is for monitoring program performance.

    This avoids the "rollover" in clockDumpable, but is less precise. */
  class secondsDumpable : 
    public fnd::dumpable<fnd::basicDumpable::dumpArg>
  {
    mzrUnit& rMzrUnit;

  public:
    secondsDumpable(mzrUnit& refMzrUnit);
    
    virtual void
    doDump(const fnd::basicDumpable::dumpArg& rDumpArg) const;
  };
}

#endif // MZR_SECONDSDUMPABLE_H

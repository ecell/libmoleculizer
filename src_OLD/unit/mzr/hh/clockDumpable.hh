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

#ifndef MZR_CLOCKDUMPABLE_H
#define MZR_CLOCKDUMPABLE_H

#include "fnd/dumpable.hh"

namespace mzr
{
  /*! \ingroup dumpGroup
    \brief Dumps dumps clock time since program start in fractions of a second.

    This is for monitoring program performance.

    The operating system reports the clock time as an integer number
    of hundredths of a second.  Eventually, this integer number
    "rolls over" and becomes negative.  This limits the usefulness
    of clockDumpable. */
  class clockDumpable :
    public fnd::dumpable<fnd::basicDumpable::dumpArg>
  {
  public:
    clockDumpable(void);
    
    virtual void
    doDump(const fnd::basicDumpable::dumpArg& rDumpArg) const;
  };
}

#endif // MZR_CLOCKDUMPABLE_H

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#include "utl/utility.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/moleculizer.hh"
#include "mzr/tabDumpEvent.hh"

namespace mzr
{
tabDumpEvent::
tabDumpEvent(double dumpPeriod,
const std::string& fileName) :
fnd::dumpStream(fileName),
period(dumpPeriod)
{}

fnd::eventResult
tabDumpEvent::happen(moleculizer& rMolzer)
throw(std::exception)
{
// Dump to output file.
doDump();

// Reschedule the event after time period.
//     rMolzer.eventQ.scheduleEvent(this,
// 				 rMolzer.eventQ.getSimTime() + period);

return fnd::go;
}

void
tabDumpEvent::
insertTaggedDumpStreamElts(xmlpp::Element* pParent) const
throw(std::exception)
{
xmlpp::Element* pTaggedDumpStreamElt
= pParent->add_child(eltName::taggedDumpStream);

pTaggedDumpStreamElt
->set_attribute(eltName::taggedDumpStream_fileNameAttr,
getFileName());

pTaggedDumpStreamElt
->set_attribute(eltName::taggedDumpStream_dumpPeriodAttr,
utl::stringify<double>(period));

std::for_each(speciesStreams.begin(),
speciesStreams.end(),
std::bind2nd
(std::mem_fun
(&mzrSpeciesStream::insertTaggedSpeciesStreamRef),
pTaggedDumpStreamElt));
}
}

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

#include <algorithm>
#include "mzr/moleculizer.hh"
#include "mzr/event.hh"
#include "mzr/species.hh"
#include "mzr/reaction.hh"
#include "mzr/mzrUnit.hh"

namespace mzr
{
  double
  event::getTime(void)
  {
    if(scheduled) return iQueueEntry->first;
    else return never;
  }

  // This kind of event is for initializing the simulation;
  // it just instantly creates the given number of the given
  // species.
  createEvent::createEvent(species* pSpecies,
			   int count,
			   mzrUnit& refMzrUnit) :
    pSpeciesToCreate(pSpecies),
    howMany(count),
    rMzrUnit(refMzrUnit)
  {}

  event::eventResult
  createEvent::doEvent(moleculizer& rMolzer,
		       mzrUnit& rMzrUnit)
  {
    std::set<reaction*> affectedReactions;
    pSpeciesToCreate->update(affectedReactions,
			     howMany);

    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     scheduleReaction(rMolzer.eventQ,
			      rMzrUnit));

    return go;
  }

  event::eventResult
  volumeEvent::doEvent(moleculizer& rMolzer,
		       mzrUnit& rMzrUnit)
  {
    std::set<reaction*> affectedReactions;
  
    // Reset the volume.
    rMzrUnit.getMolarFactor().updateVolume(vol,
					   affectedReactions);

    // Reschedule everybody.  This invocation stinks to high
    // heaven of confused reference structure in this program.
    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     scheduleReaction(rMolzer.eventQ,
			      rMzrUnit));
    return go;
  }

  event::eventResult
  growEvent::doEvent(moleculizer& rMolzer,
		     mzrUnit& rMzrUnit)
  {

    double currentVolume = rMzrUnit.getMolarFactor().getVolume();
  
    std::set<reaction*> affectedReactions;
    rMzrUnit.getMolarFactor().updateVolume(currentVolume * factor,
					   affectedReactions);

    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     scheduleReaction(rMolzer.eventQ,
			      rMzrUnit));

    // Reschedule this event for later.
    eventQueue& rQ = rMolzer.eventQ;
    rQ.scheduleEvent(this,
		     rQ.getSimTime() + period);

    return go;
  }

  event::eventResult
  stopEvent::doEvent(moleculizer& rMolzer,
		     mzrUnit& rMzrUnit)
  {
    return stop;
  }

  tabDumpEvent::tabDumpEvent(double dumpPeriod,
			     const std::string& fileNameString) :
    pOs(0),
    pFileStream(0),
    fileName(fileNameString),
    period(dumpPeriod)
  {}

  // Initializes output stream and writes column headers.
  void
  tabDumpEvent::init(void)
  {
    if(fileName == "-") pOs = &std::cout;
    else if(fileName == "+") pOs = &std::cerr;
    else
      {
	// Open the file for appending, so that if a simulation is continued,
	// the output will also be continued instead of starting over at the
	// second start time.
	pFileStream = new std::ofstream(fileName.c_str(), std::ios_base::app);

	if(! (*pFileStream))
	  throw badDumpFileXcpt(fileName);

	pOs = pFileStream;
      }

    // The header line is a "comment line" for gnuplot.
    (*pOs) << "#";
    
    // Emit the headers as a tab-separated list.
    std::vector<dumpable*>::const_iterator ipDumpable = dumpables.begin();
    if(ipDumpable != dumpables.end())
      {
	(*ipDumpable)->dumpHeader(*pOs);

	while(++ipDumpable != dumpables.end())
	  {
	    (*pOs) << '\t';
	    (*ipDumpable)->dumpHeader(*pOs);
	  }
      }

    // Emit a newline.
    (*pOs) << std::endl;
  }

  void
  tabDumpEvent::insertTaggedDumpStreamElts(xmlpp::Element* pParent) const
    throw(std::exception)
  {
    xmlpp::Element* pTaggedDumpStreamElt
      = pParent->add_child(eltName::taggedDumpStream);

    pTaggedDumpStreamElt
      ->set_attribute(eltName::taggedDumpStream_fileNameAttr,
		      fileName);

    pTaggedDumpStreamElt
      ->set_attribute(eltName::taggedDumpStream_dumpPeriodAttr,
		      domUtils::stringify<double>(period));

    std::for_each(speciesDumpables.begin(),
		  speciesDumpables.end(),
		  std::bind2nd
		  (std::mem_fun
		   (&speciesDumpable::insertTaggedSpeciesStreamRef),
		   pTaggedDumpStreamElt));
  }

  event::eventResult
  tabDumpEvent::doEvent(moleculizer& rMolzer,
			mzrUnit& rMzrUnit)
  {
    // Emit the dumpables as a tab-separated list.  Not doing this
    // functionally because the last dumpable is a special case.
    std::vector<dumpable*>::iterator ipDumpable = dumpables.begin();
    while(ipDumpable != dumpables.end())
      {
	// Emit a dumpable.
	(*ipDumpable)->doDump(*pOs);

	// Emit a tab unless we're at the end of the line.
	if(++ipDumpable != dumpables.end())
	  {
	    (*pOs) << '\t';
	  }
      }

    // Emit a newline.
    (*pOs) << std::endl;
      
    // Reschedule the event after time period.
    rMolzer.eventQ.scheduleEvent(this,
				 rMolzer.eventQ.getSimTime() + period);

    return go;
  }

  // Dumps state to a file for export to other programs or other processing.
  event::eventResult
  dumpStateEvent::doEvent(moleculizer& rMolzer,
			  mzrUnit& rMzrUnit) throw(std::exception)
  {
    xmlpp::Document* pOutputDoc = rMolzer.makeDomOutput();

    pOutputDoc->write_to_file(outputFileName);

    delete pOutputDoc;

    return go;
  }


  // Stops reaction generation.
  event::eventResult
  noReactEvent::doEvent(moleculizer& rMolzer,
			mzrUnit& rMzrUnit)
  {
    rMzrUnit.generateOk = false;
    return go;
  }
}

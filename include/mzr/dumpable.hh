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

#ifndef DUMPABLE_H
#define DUMPABLE_H

#include "domUtils/domUtils.hh"
#include "mzr/mzrEltName.hh"

namespace mzr
{
  class species;
  class mzrUnit;
  class moleculizer;

  /*! \defgroup dumpGroup Dumping and tracing
    \ingroup mzrGroup
    \brief Dumping populations and tracing reactions. */

  /*! \file dumpable.hh
    \ingroup dumpGroup
    \brief Defines dumpable base class and built-in descendents.

    This file gives the abstract dumpable base class, the interface
    for dumping output when a dumpEvent happens.

    The most common simulation data to dump this way are species
    populations; each species is enabled to dump in this way.  Most of
    the descendents of dumpable defined in this file are for dumping
    other simulation data. */

  /*! \ingroup dumpGroup
    \brief Abstract base class for dumpables. */
  class dumpable
  {
    // Used as a column header in the dump file.
    std::string name;
    
  public:
    dumpable(const std::string& rName) :
      name(rName)
    {}
    
    virtual
    ~dumpable(void)
    {}

    const std::string&
    getName(void) const
    {
      return name;
    }

    virtual void
    doDump(std::ostream& rOs) const = 0;

    void
    dumpHeader(std::ostream& rOs) const
    {
      rOs << getName();
    }
  };

  // This class is for dumpables that actually dump the population of one or
  // more species.  It exists so that it can list the species whose populations
  // it dumps for export to other programs.
  class speciesDumpable : public dumpable
  {
  public:
    speciesDumpable(const std::string& rName) :
      dumpable(rName)
    {}
    
    // For state dump.

    virtual void
    insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
      throw(std::exception) = 0;

    void
    insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
      throw(std::exception);
  };

  /*! \ingroup dumpGroup
    \brief Dumps the simulation time.

    Automatically used as the first field by most dump commands. */
  class simTimeDumpable : public dumpable
  {
    mzr::moleculizer& rMolzer;
    
  public:
    simTimeDumpable(mzr::moleculizer& rMoleculizer) :
      dumpable(eltName::statStream_simTime),
      rMolzer(rMoleculizer)
    {}
      
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps the population of a species. */
  class singleSpeciesDumpable : public speciesDumpable
  {
    const species& rSpecies;
  public:
    singleSpeciesDumpable(const std::string& rName,
			  const species& rSpeciesToDump) :
      speciesDumpable(rName),
      rSpecies(rSpeciesToDump)
    {}

    virtual void
    doDump(std::ostream& rOs) const;

    virtual void
    insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
      throw(std::exception);
  };

  // For species streams that will only be notified of species that they
  // should dump.  Doesn't do any filtering of the species if which it's
  // notified.
  //
  // This would be useful for dumping all the species that display a certain
  // feature; e.g. omniPlexDumpable.
  class multiSpeciesDumpable : public speciesDumpable
  {
    // Function class to sum populations of dumped species.
    class accumPop;

    // Function class to insert tag of one of the dumped species.
    class insertTag;
    
    std::set<const species*> dumpedSpecies;
  public:
    multiSpeciesDumpable(const std::string& rName) :
      speciesDumpable(rName)
    {}

    void
    addSpecies(const species* pSpecies)
    {
      dumpedSpecies.insert(pSpecies);
    }

    virtual void
    doDump(std::ostream& rOs) const;

    virtual void
    insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
      throw(std::exception);
  };

  /*! \ingroup dumpGroup
    \brief Dumps dumps clock time since program start in fractions of a second.

    This is for monitoring program performance.

    The operating system reports the clock time as an integer number
    of hundredths of a second.  Eventually, this integer number
    "rolls over" and becomes negative.  This limits the usefulness
    of clockDumpable. */
  class clockDumpable : public dumpable
  {
  public:
    clockDumpable(void) :
      dumpable(eltName::statStream_clockTime)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps dumps clock time since program start in whole seconds.

    This is for monitoring program performance.

    This avoids the "rollover" in clockDumpable, but is less precise. */
  class secondsDumpable : public dumpable
  {
    mzrUnit& rMzrUnit;
  public:
    secondsDumpable(mzrUnit& refMzrUnit) :
      dumpable(eltName::statStream_clockSeconds),
      rMzrUnit(refMzrUnit)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps the number of reactions performed so far.

    This is for monitoring program performance and progress.

    Reaction complexity varies considerably, so that number of reactions
    performed may not correlate well with amount of work done.

    \sa activationCountDumpable
  */
  class reactEventCountDumpable : public dumpable
  {
  public:
    reactEventCountDumpable(void) :
      dumpable(eltName::statStream_reactionEventCount)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps the number of reactions that have occurred at least once.
  */
  class reactCountDumpable : public dumpable
  {
  public:
    reactCountDumpable(void) :
      dumpable(eltName::statStream_reactionCount)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps the number of reaction activations so far.

    "Activation" means propensity calculation and rescheduling.  One
    reaction event may cause recalculation of propensity and
    rescheduling of many or few other reactions.  This recalculation and
    rescheduling is the bulk of the work, so that the amount of work in
    one reaction may vary considerably.  The activation count should be
    a better indicator of the amount of work done by the simulation so
    far. */
  class activationCountDumpable : public dumpable
  {
  public:
    activationCountDumpable(void) :
      dumpable(eltName::statStream_activationCount)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps the current simulation volume. */
  class volumeDumpable : public dumpable
  {
    mzrUnit& rMzrUnit;
    
  public:
    volumeDumpable(mzrUnit& refMzrUnit) :
      dumpable(eltName::statStream_volume),
      rMzrUnit(refMzrUnit)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };

  /*! \ingroup dumpGroup
    \brief Dumps the number of species that have been encountered. */
  class speciesCountDumpable : public dumpable
  {
  public:
    speciesCountDumpable(void) :
      dumpable(eltName::statStream_speciesCount)
    {}
    
    virtual void
    doDump(std::ostream& rOs) const;
  };
}

#endif

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

#ifndef BASICSPECIES_H
#define BASICSPECIES_H

#include <sstream>
#include "fnd/physConst.hh"
#include "utl/xcpt.hh"
#include "utl/utility.hh"
#include "fnd/stateVar.hh"
#include "fnd/notifier.hh"

namespace fnd
{
    template<class reactionType>
    class basicSpecies :
                public stateVar<int, reactionType>,
                public onceNotifier
    {
    class negativePopulationXcpt :
                    public utl::xcpt
        {
            static std::string
            mkMsg (basicSpecies* pBasicSpecies,
                   int delta,
                   int newValue)
            {
                typename std::ostringstream msgStream;
                msgStream << "Species "
                << pBasicSpecies
                << " assumed negative value "
                << newValue
                << " after update by delta "
                << delta
                << ".";
                return msgStream.str();
            }
        public:
            negativePopulationXcpt (basicSpecies* pBasicSpecies,
                                    int delta,
                                    int newValue) :
                    utl::xcpt (mkMsg (pBasicSpecies,
                                      delta,
                                      newValue) )
            {}
        };

        static int speciesCount;

    public:

        basicSpecies (int initialPop = 0) :
                stateVar<int, reactionType> (initialPop)
        {
            speciesCount++;
        }

        virtual
        ~basicSpecies (void)
        {}

        static int
        getSpeciesCount (void)
        {
            return speciesCount;
        }

// Gives this address as a hex string.
        typename std::string
        getTag (void) const
        {
            return utl::stringify<const basicSpecies*> (this);
        }

// For possibly getting a more humanly-readable, informative name.
        virtual typename std::string
        getName (void) const
        {
            return getTag();
        }

        int
        getPop (void) const
        {
            return this->getValue();
        }

// Note that this should only be used for initialization, since
// it doesn't get help keep track of the reactions whose propensities
// will change because of the change in population of this species.
        void
        setPop (int newPop)
        {
            this->setValue (newPop);
        }

        double
        getConcentration (double volume) const
        {
            return ( (double) getPop() ) / (fnd::avogadrosNumber * volume);
        }

// For updating the population using a positive or negative
// delta.  Throws an exception if the population goes negative.
        void
        update (int delta,
                sensitivityList<reactionType>& rAffectedReactions,
                int notificationDepth)
        throw (utl::xcpt)
        {
// This has to come first, so that the sensitivity list
// is correct for the update.
            ensureNotified (notificationDepth);

            int newValue = this->getValue() + delta;

            if (newValue < 0) throw negativePopulationXcpt (this,
                        delta,
                        newValue);
            updateValue (newValue,
                         rAffectedReactions);
        }
    };

    template<class reactionType>
    int basicSpecies<reactionType>::speciesCount = 0;
}

#endif // BASICSPECIES_H

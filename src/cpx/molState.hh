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
// along with Moleculizer; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001

//                     Email: lok@molsci.org
//   
/////////////////////////////////////////////////////////////////////////////

#ifndef MOLSTATE_H
#define MOLSTATE_H

#include "utl/debug.hh"

namespace cpx
{
    /*! \ingroup plexSpeciesGroup
      \brief (Concrete) base class for state of molecules. */
    class molState
    {
    protected:
        // Other contributions from molecular state may go into the true
        // molecular weight; hence the name.
        double baseWeight;

    public:
        molState(double molWeight) :
            baseWeight(molWeight)
        {}
        virtual ~molState()
        {}
        virtual double
        getMolWeight(void) const
        {
            return baseWeight;
        }

        bool
        operator<(const molState& rRight) const
        {
            return baseWeight < rRight.baseWeight;
        }

        // DEBUG
        virtual 
        string
        getName() const
        {
            ostringstream oss;
            oss << "mol state with weight '" 
                << baseWeight 
                << "'";
            return oss.str();

        }
        

    };

    // Now not sure if this was such a good idea...
    typedef const molState* molParam;
}

#endif // MOLSTATE_H

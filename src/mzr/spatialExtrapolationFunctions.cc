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
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2008
//
// Modifing Authors:
//
//

#include "spatialExtrapolationFunctions.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrReaction.hh"
#include <cmath>

namespace mzr
{

    double getProteinDiffusionCoef()
    {
        // Units returned are in micrometers^2/sec.
        return proteinDiffusionCoeff;
    }

    double getSmallMolDiffusionCoef()
    {
        // Units returned are in micrometers^2/sec.
        return smallMolDiffusionCoeff;
    }


    void setProteinDiffusionCoef(double rate)
    {
        // Units must be in micrometer^2/sec.
        proteinDiffusionCoeff = rate;
    }

    void setSmallMolDiffusionCoef(double rate)
    {
        // Units must be in micrometer^2/sec.
        smallMolDiffusionCoeff = rate;
    }

    double calculateSumOfRadii( const mzr::mzrReaction* pRxn )
    {
        double theSum = 0.0f;

        typedef std::map<mzr::mzrSpecies*, int> multMap;
        BOOST_FOREACH(const multMap::value_type& vt, pRxn->getReactants() )
        {
            // theSum += multiplicity * 
            theSum += vt.second * extrapolateMolecularRadius( vt.first );
        }

        return theSum;
    }

    double getDiffusionCoeffFromSpecies( const mzr::mzrSpecies* pSpecies)
    {
        // This is kind of hacky but probably works.  I am using 1500 daltons
        // as the cutoff for small molecules and proteins.  

        if ( pSpecies->getWeight() > 1500 )
        {
            return getSmallMolDiffusionCoef();
        }
        else
        {
            return getProteinDiffusionCoef();
        }
    }

    double calculateSumOfDiffusionCoefficients( const mzr::mzrReaction* pRxn )
    {
        double theSum = 0.0f;

        typedef std::map<mzr::mzrSpecies*, int> multMap;
        BOOST_FOREACH(const multMap::value_type& vt, pRxn->getReactants() )
        {
            // theSum += multiplicity * particle-type specific diffusion coeff.
            theSum += vt.second * getDiffusionCoeffFromSpecies( vt.first );
        }

        return theSum;
    }

    double extrapolateIntrinsicReactionRate(const mzr::mzrReaction* pRxn)
    {
        static const double FourPi = 4.0f * 3.141592653859;

        // We have to solve the equation 1/k = 1/kA + 1/ kD for kA, where 
        double k = pRxn->getRate();
        double kD = FourPi * calculateSumOfRadii(pRxn) * calculateSumOfDiffusionCoefficients( pRxn );
        
        return ( k * kD ) / ( kD - k );
    }
 
    double extrapolateMolecularRadius( const mzr::mzrSpecies* pSpecies)
    {
        // This number is obtained by converting to angstroms to m^3 by using the conversion
        // factor 1.22 g/cm^3 and then solving for R in the equation V = 4/3 * Pi * R^3
        const double Factor( 3.0774371428132849e30 );
        const double radical( 1.0 / 3.0 );

        double newSpeciesConvertedMass = pSpecies->getWeight() / Factor;
        return std::pow( newSpeciesConvertedMass, radical );
    }


    double extrapolateMolecularRadius(const double& mass)
    {
        // This number is obtained by converting to angstroms to m^3 by using the conversion
        // factor 1.22 g/cm^3 and then solving for R in the equation V = 4/3 * Pi * R^3
        const double Factor( 3.0774371428132849e30 );
        const double radical( 1.0 / 3.0 );

        double newSpeciesConvertedMass = mass / Factor;
        return std::pow( newSpeciesConvertedMass, radical );

    }

}


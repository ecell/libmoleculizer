///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


#include "moleculizer.hh"

namespace mzr
{

void
moleculizer::installKForNewUnaryReaction( const mzrReaction* mzrReaction )
{
    // Get the originating reaction generator.  Look up its reaction rate in the
    // appropriate moleculizer chart and use that to install the new K for this reaction
    // in the appropriate chart.

    if ( mzrReaction->getNumberOfReactants() != 1 ) return;

    const fnd::coreRxnGen* parentReactionGen( mzrReaction->getOriginatingRxnGen() );
    double basicUnaryReactionRate = unaryReactionRatesParameterLookup[ parentReactionGen ];

    reactionRateChart.insert( std::make_pair( mzrReaction, basicUnaryReactionRate ) );
}

void
moleculizer::installKForNewBinaryReaction( const mzrReaction* mzrReaction )
{
    if ( mzrReaction->getNumberOfReactants() != 2 ) return;

    const fnd::coreRxnGen* parentReactionGen( mzrReaction->getOriginatingRxnGen() );

    double basicForwardReactionRate = binaryReactionRatesParameterLookup[ parentReactionGen ];
    reactionRateChart.insert( std::make_pair( mzrReaction, basicForwardReactionRate ) );
}

void
moleculizer::installKaForNewBinaryReaction( const mzrReaction* mzrReaction )
{
    if ( mzrReaction->getNumberOfReactants() != 2 ) return;

    const fnd::coreRxnGen* parentGen( mzrReaction-> getOriginatingRxnGen() );

    double basicForwardActivationEnergy = binaryActivationEnergiesParameterLookup[parentGen];

    // We do nothing to extrapolate ( I don't think there is anything *to* do with it)
    // and just install it in the right place.

    binaryActivationEnergyChart.insert( std::make_pair( mzrReaction, basicForwardActivationEnergy ) );

}

void
moleculizer::installRadiusForNewSpecies( const mzrReaction* mzrReaction )
{
    // For each of the new species (products in this new reaction), calculate its
    // radius, assuming the protein is globular, and based on a protein density of
    // 1.22g/cm^3.
    // The value type here is a pair(mzrSpecies*, int) -- int represents the multiplicity.

    BOOST_FOREACH( const mzrReaction::multMap::value_type& vt, mzrReaction->getProducts() )
    {
        // I'm not certain, but I'm reasonably sure this can happen, for instance in
        // the decomposition of a species into a known and an unknown product.
        if ( radiusChart.find( vt.first->getName() ) != radiusChart.end() )
        {
            ; // Do nothing.
        }
        else
        {
            std::string speciesName( vt.first->getName() );
            double newRadius = calculateNewRadiusForSpecies( vt.first );

            std::cout << "New species " << speciesName << " radius=" << newRadius << "m" << std::endl;
            this->radiusChart.insert( std::make_pair( speciesName, newRadius ) );
        }
    }

}

double
moleculizer::calculateNewRadiusForSpecies( const mzrSpecies* ptrSpecies )
{
    // This number is obtained by converting to angstroms to m^3 by using the conversion
    // factor 1.22 g/cm^3 and then solving for R in the equation V = 4/3 * Pi * R^3
    const double Factor( 3.0774371428132849e30 );
    const double radical( 1.0 / 3.0 );

    double newSpeciesConvertedMass = ptrSpecies->getWeight() / Factor;
    return std::pow( newSpeciesConvertedMass, radical );
}

void
moleculizer::installKdForNewSpecies( const mzrReaction* mzrReaction )
{
    // The diffusion coefficient for all new species should be (for now)
    // the min/average/max of all the product species.

    std::vector<double> kdOfReactants;

    BOOST_FOREACH( const mzrReaction::multMap::value_type& vt, mzrReaction->getReactants() )
    {

        double reactantKd;
        if ( k_DChart.find( vt.first->getName() )== k_DChart.end() )
        {
            throw missingExtrapolationParameter( vt.first->getName(),
                                                 "k_D" );
        }
        else
        {
            reactantKd = k_DChart.find( vt.first->getName() )->second;
        }

        std::cout << "Reactant " << vt.first->getName() << " has kD " << reactantKd << std::endl;
        kdOfReactants.push_back( k_DChart[ vt.first->getName()] );
    }

    std::cout << "New Species being created: " << kdOfReactants.size() << std::endl;

    double minimumVal = *std::min_element( kdOfReactants.begin(),
                                           kdOfReactants.end() );

    std::cout << "Product has kd=" << minimumVal << std::endl;


    BOOST_FOREACH( const mzrReaction::multMap::value_type& vt, mzrReaction->getProducts() )
    {
        if ( k_DChart.find( vt.first->getName() ) != k_DChart.end() ) continue;

        std::cout << "New species " << vt.first->getName() << " kD=" << minimumVal << "m" << std::endl;

        k_DChart[ vt.first->getName()] = minimumVal;
    }


}

}




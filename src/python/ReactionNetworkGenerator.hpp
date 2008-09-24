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
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#ifndef REACTIONNETWORKGENERATOR_HPP
#define REACTIONNETWORKGENERATOR_HPP

#include "mzr/moleculizer.hh"
#include "mzr/mzrException.hh"
#include "fnd/basicReaction.hh"
#include "mzr/mzrSpecies.hh"
#include "utl/utility.hh"
#include "nmr/complexSpecies.hh"
#include "mol/mzrMol.hh"
#include <vector>
#include <string>

class Species
{
public:
//     Species()
//         :
//         name(""),
//         mass(0.0f){}

Species(const mzr::mzrSpecies& aMzrSpecies)
:
name( aMzrSpecies.getName() ),
mass( aMzrSpecies.getWeight() )
{}

std::string getName() const { return name;}
float getMass() const {return mass;}

private:
std::string name;
float mass;
};

class Reaction
{
public:
typedef fnd::basicReaction<mzr::mzrSpecies> CoreRxnType;

//    Reaction(){}
Reaction(const CoreRxnType& aReaction);

float getRate() const
{
return rate;
}

std::vector<Species>
getSubstrates() const
{
return substrates;
}

std::vector<Species>
getProducts() const
{
return products;
}

protected:
float rate;
std::vector<Species> substrates;
std::vector<Species> products;
};


struct BasicComplexRepresentation
{
typedef std::pair<int, std::string> HalfBinding;
typedef std::pair<HalfBinding, HalfBinding> BindingType;
typedef std::pair< int, std::pair<std::string, std::string> > ModType;

void
addMolNameToComplex(const std::string& molName);

void addModificationToComplex(int molIndex,
const std::string& modificationSiteName,
const std::string& modificationValue);

void addBindingToComplex( int firstMolNdx,
const std::string& bindingName1,
int secondMolNdx,
const std::string& bindingName2);

std::vector<std::string> mols;
std::vector<BindingType > bindings;
std::vector<ModType> modifications;

};


class ReactionNetworkGenerator
{

public:
ReactionNetworkGenerator()
{
ptrMoleculizer = new mzr::moleculizer;
}

~ReactionNetworkGenerator()
{
delete ptrMoleculizer;
}

void
runInteractiveMode()
{
ptrMoleculizer->RunInteractiveDebugMode();
return;
}

int
getNumUnary() const
{
ptrMoleculizer->unaryReactionList.size();
}

int getNumBinary() const
{
ptrMoleculizer->binaryReactionList.size();
}

void showDeadSpecies() const
{
ptrMoleculizer->DEBUG_showDeadSpecies();
}

void showLiveSpecies() const
{
ptrMoleculizer->DEBUG_showLiveSpecies();
}

void addRules(const std::string& filename) throw( mzr::BadRulesDefinitionXcpt )
{
try
{
delete ptrMoleculizer;
ptrMoleculizer = new mzr::moleculizer;
ptrMoleculizer->attachFileName( filename );
}
catch(mzr::BadRulesDefinitionXcpt e)
{
e.warn();

if(ptrMoleculizer)
{
delete ptrMoleculizer;
ptrMoleculizer = NULL;
}

throw e;
}
catch(utl::xcpt x)
{
delete ptrMoleculizer;
x.warn();
throw mzr::BadRulesDefinitionXcpt();
}
catch(...)
{
delete ptrMoleculizer;
utl::xcpt x("Unknown Error while adding rules. (FDLLAE)");
x.wailAndBail();
}
}

std::vector<Reaction>
getBinaryReactions(const std::string& species1,
const std::string& species2) throw( mzr::IllegalNameXcpt );

std::vector<Reaction>
getUnaryReactions(const std::string& species1) throw (mzr::IllegalNameXcpt );

Species
getSpecies(const std::string& species)throw (mzr::IllegalNameXcpt );

bool
checkSpeciesNameLegality( const std::string& species1) throw (mzr::IllegalNameXcpt );


// For now I am writing two versions of this function.  This 'stricter' one insists on
// molecules in the ComplexRepresentation being a priori defined in the moleculizer
// rules.  This is so that "default" modifications that aren't in the complex representation.
// But that might just be stupid.
std::string
generateNameFromBasicComplexRepresentation(const BasicComplexRepresentation& aBCR);

std::string
generateNameFromBasicComplexRepresentationStrict(const BasicComplexRepresentation& aBCR);

void showAllReactions();
void incrementSpecies(const std::string& species);
int getNumberOfSpecies();
int getNumberOfReactions();

void
showAllSpecies()
{
ptrMoleculizer->DEBUG_showAllSpecies();
}


protected:
mzr::moleculizer* ptrMoleculizer;
};

class modificationNotOfNdx
{
public:
int ndx;
modificationNotOfNdx( int index)
:
ndx( index )
{}

bool operator()(const BasicComplexRepresentation::ModType& theModType)
{
return theModType.first != ndx;
}
};


#endif

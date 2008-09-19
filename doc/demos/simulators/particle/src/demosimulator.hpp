/////////////////////////////////////////////////////////////////////////////
// libMoleculizer - a species and reaction generator for reaction networks
// Copyright (C) 2008  Molecular Sciences Institute
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
// This file was authored by Nathan Addy <addy@molsci.org> 
// Contact information:
//   Nathan Addy, Scientific Programmer     Email: addy@molsci.org
//   The Molecular Sciences Institute
//   2168 Shattuck Ave.
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////


#include "mzr/moleculizer.hh"
#include <string>
#include <map>

class SimplifiedStochasticTypeSimulator
{
public:
    SimplifiedStochasticTypeSimulator(std::string modelfile, std::string rulesfile);

    void singleStep();
    void attachRuleFile(std::string rulesfile);
    bool assertModelValidity();
    void engageModel();


    mzr::moleculizer speciesReactionGenerator;
    std::map<std::string, int> theModel;
    std::vector<mzr::mzrReaction> theReactionList;

public:
    typedef std::pair<std::string, int> modelPairType;

private:
    void createModelFromFile(const std::string& modelFile, std::map<std::string, int>& model);
    
    
        
        
};





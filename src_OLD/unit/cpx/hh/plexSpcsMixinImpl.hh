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

#ifndef CPX_PLEXSPECIESIMPL_H
#define CPX_PLEXSPECIESIMPL_H


namespace cpx
{
  class addMolWeight :
    public std::unary_function<molParam, void>
  {
    double& rTotal;
  public:
    addMolWeight(double& refTotalMolWeight) :
      rTotal(refTotalMolWeight)
    {
    }

    void
    operator()(molParam pState) const
    {
      rTotal += pState->getMolWeight();
    }
  };
  
  template<class plexFamilyT>
  double
  plexSpeciesMixin<plexFamilyT>::
  getWeight(void) const
  {
    double totalMolWeight = 0.0;
    
    std::for_each(molParams.begin(),
		  molParams.end(),
		  addMolWeight(totalMolWeight));

    return totalMolWeight;
  }
  
  template<class plexFamilyT>
  std::string
  plexSpeciesMixin<plexFamilyT>::
  getInformativeName(void) const
  {
    // Let's start with the names of the mols in the order that
    // they appear in the pardigm, separated by colons.
    //
    // It might be nice to traverse the mols in some kind of "connectivity
    // order" later on, if called for.

    // Bad that the type "gets out" like this; could be avoided, of course,
    // by using lots more typedefs.
    typedef typename std::vector<typename plexFamilyT::molType*> molVectorType;
    
    const molVectorType& rMols
      = rFamily.getParadigm().mols;

    std::string theName;

    typename molVectorType::const_iterator iMol = rMols.begin();
    if(rMols.end() != iMol)
      {
	typename plexFamilyT::molType* pMol = *iMol++;
	theName = pMol->getName();
      }

    while(rMols.end() != iMol)
      {
	typename plexFamilyT::molType* pMol = *iMol++;
	theName += "_";
	theName += pMol->getName();
      }

    return theName;
  }
}

#endif // CPX_PLEXSPECIESIMPL_H

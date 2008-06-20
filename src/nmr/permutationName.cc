#include "permutationName.hh"
#include <iostream>
using std::cout;
using std::endl;

namespace nmr
{
    bool PermutationName::operator<(PermutationNameCref aPermutationName) const
    {
        counter++;
//         if (counter==1485038)
//         {
//             std::cout << "PTL1:\t" << theCorrespondingPartialTokenList << std::endl;
//             // std::cout << "PTL2:\t" << aPermutationName.theCorrespondingPartialTokenList << std::endl;
//         }
        
        return (theCorrespondingPartialTokenList < aPermutationName.theCorrespondingPartialTokenList);
    }

    void PermutationName::print() const
    {
        std::cout << "[";
        for (unsigned int ndx = 0; ndx != thePermutation.getDimension(); ++ndx)
        {
            std::cout << thePermutation[ndx];
            if (ndx != thePermutation.getDimension() - 1)
            {
                std::cout << ", ";
            }

        }
        cout << "]" << endl;
    }
        
    int PermutationName::counter = 0;
}


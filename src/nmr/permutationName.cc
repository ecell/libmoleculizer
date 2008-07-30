#include "permutationName.hh"

namespace nmr
{
    bool PermutationName::operator<(PermutationNameCref aPermutationName) const
    {
        return (theCorrespondingPartialTokenList < aPermutationName.theCorrespondingPartialTokenList);
    }

    int PermutationName::counter = 0;
}


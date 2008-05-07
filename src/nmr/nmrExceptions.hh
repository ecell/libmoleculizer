/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008  Nathan Addy
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
//   Nathan Addy, Research Associate     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef NMR_NMREXCEPTIONS_HH
#define NMR_NMREXCEPTIONS_HH

#include "utl/xcpt.hh"

namespace nmr
{

    class GeneralNmrXcpt
        :
        public utl::xcpt
    {
        public:
        GeneralNmrXcpt(const std::string& message)
            :
            utl::xcpt( message )
        {}
    };

    class MissingMolAliasXcpt
        :
        public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& unfoundMolAlias)
        {
            std::ostringstream oss;
            oss << "Error: mol alias '"
                << unfoundMolAlias 
                << "' was looked up but not found in complex.";
            return oss.str();
        }

    public:
        MissingMolAliasXcpt(const std::string& unfoundMolAlias) 
            : 
            utl::xcpt( mkMsg( unfoundMolAlias) )
        {}
    };


    class MissingBindingSiteXcpt
        :
        public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& missingBindingSite, const std::string& ostensiblyOwningMol)
        {
            std::ostringstream oss;
            oss << "Error: Binding site '"
                << missingBindingSite
                << "' was looked up on mol '"
                << ostensiblyOwningMol 
                << "' but was not found.";
            return oss.str();
        }

    public:
        MissingBindingSiteXcpt(const std::string& missingBindingSite, 
                            const std::string& ostensiblyOwningMol) 
            : 
            utl::xcpt( mkMsg( missingBindingSite, ostensiblyOwningMol) )
        {}
    };

    class DuplicateMolAliasXcpt
        :
        public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& molType, const std::string& molAlias )
        {
            std::ostringstream oss;
            oss << "Error: mol of type '"
                << molType 
                << "' was attempted to be added to a complex under a duplicate alias '"
                << molAlias
                << "'";
            return oss.str();
        }

    public:
        DuplicateMolAliasXcpt( const std::string& molType, const std::string& molAlias )
            :
            utl::xcpt( mkMsg(molType, molAlias) )
        {}
    };
    

    class BadPermutationConstructorXcpt
        :
        public utl::xcpt
    {
        static std::string
        mkMsg( int ActualValue, unsigned int pos)
        {
            std::ostringstream oss;
            oss << "Error in permutation constructor.  The permutation being extended should "
                << "have UNDEFINED value at position = '" 
                << pos
                << "'; however value("
                << pos 
                << ") = " 
                << ActualValue << ".";
            return oss.str();
        }
    public:

        BadPermutationConstructorXcpt( int ActualValue, unsigned int pos)
            :
            utl::xcpt(mkMsg(ActualValue, pos) )
        {}
    };
    
    class IncompatiblePermutationsXcpt
        :
        public utl::xcpt
    {
        static std::string mkMsg( unsigned int dimension1, unsigned int dimension2)
        {
            std::ostringstream oss;
            oss << "Error: Two permutations with incompatible dimensionalities ('" 
                << dimension1
                << "' and '"
                << dimension2
                << "') were used together.";
            return oss.str();
        }

    public:
        IncompatiblePermutationsXcpt( unsigned int dimension1, unsigned int dimension2)
            :
            utl::xcpt( mkMsg( dimension1, dimension2) )
        {}
    };

    class BadPermutationIndexXcpt
        :
        public utl::xcpt
    {
        static std::string
        mkMsg( unsigned int dimension, unsigned int badNdx)
        {
            std::ostringstream oss;
            oss << "Bad domain element used with permutation.  Permutation has dimension " 
                << dimension
                << ", however the value at "
                << badNdx 
                << " was called for.";

            return oss.str();
        }

    public:
        BadPermutationIndexXcpt( unsigned int permutationDimension, unsigned int badNdx)
            :
            utl::xcpt( mkMsg(permutationDimension, badNdx) )
        {}
    };

    class NoSuchNameEncoderXcpt 
        :
        public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& rBadNameEncoderName)
        {
            std::ostringstream msgStream;
            msgStream << "There is no NameEncoder provided in this distribution with the name '"
                      << rBadNameEncoderName 
                      << "'.";
            return msgStream.str();
        }
    
    public:
        NoSuchNameEncoderXcpt( const std::string& rBadName):
            utl::xcpt( mkMsg(rBadName) )
        {}
    };
    class NoSuchBindingSiteXcpt : public utl::xcpt
    {
        static std::string 
        mkMsg( const std::string& rNonexistentBindingSiteName )
        {
            std::ostringstream msgStream;
            msgStream << "BindingSite '" 
                      << rNonexistentBindingSiteName 
                      << "' does not exist.";
            return msgStream.str();
        }

        static std::string 
        mkMsg( const std::string& rMolName, const std::string& rNonexistentBindingSiteName )
        {
            std::ostringstream msgStream;
            msgStream << "BindingSite '" 
                      << rNonexistentBindingSiteName 
                      << "' does not exist on mol '"
                      << rMolName
                      << "'";

            return msgStream.str();
        }

    public:
        NoSuchBindingSiteXcpt(const std::string& rNonexistentBindingSiteName)
            :
            utl::xcpt( mkMsg(rNonexistentBindingSiteName) )
        {}

        NoSuchBindingSiteXcpt(const std::string& rMolName, const std::string& rNonexistentBindingSiteName)
            :
            utl::xcpt( mkMsg( rMolName, rNonexistentBindingSiteName) )
        {}
    };

    class NoSuchModificationSiteXcpt : public utl::xcpt
    {
        static std::string 
        mkMsg( const std::string& rNonexistentModificationSiteName )
        {
            std::ostringstream msgStream;
            msgStream << "ModificationSite '" 
                      << rNonexistentModificationSiteName 
                      << "' does not exist.";
            return msgStream.str();
        }

        static std::string 
        mkMsg( const std::string& rMolName, const std::string& rNonexistentModificationSiteName )
        {
            std::ostringstream msgStream;
            msgStream << "ModificationSite '" 
                      << rNonexistentModificationSiteName 
                      << "' does not exist on mol '"
                      << rMolName 
                      << "'.";
            return msgStream.str();
        }

    public:
        NoSuchModificationSiteXcpt(const std::string& rNonexistentModificationSiteName)
            :
            utl::xcpt(mkMsg( rNonexistentModificationSiteName))
        {}

        NoSuchModificationSiteXcpt(const std::string& rMolName, const std::string& rNonexistentModificationSiteName)
            :
            utl::xcpt(mkMsg( rMolName, rNonexistentModificationSiteName))
        {}
    };
    
    class badBindingNameXcpt : public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& bindingList)
        {
            std::ostringstream msgStream;
            msgStream << "Binding phrase '" 
                      << bindingList
                      << "' could not be decoded as it is an illegal binding list.";

            return msgStream.str();
        }

    public:
        badBindingNameXcpt(const std::string& bindingName)
            :
            utl::xcpt( mkMsg(bindingName) )
        {}
    };

    class badModificationNameXcpt : public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& modificationString)
        {
            std::ostringstream msgStream;
            msgStream << "Modification string '"
                      << modificationString 
                      << "' could not be decoded as it is an illegal modification string.";
            return msgStream.str();
        }

    public:
        badModificationNameXcpt(const std::string& modificationString )
            :
            utl::xcpt( mkMsg(modificationString) )
        {}
    };


    class encodeDecodeInconsistencyXcpt : public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& encodingScheme)
        {
            std::ostringstream msgStream;
            msgStream << "ComplexSpecies encoding scheme '"
                      << encodingScheme
                      << "' does not have the property that decode( encode( \"A Complex Species State\") ) == \"A Complex Species State\"";
            return msgStream.str();
        }

    public:
        encodeDecodeInconsistencyXcpt(const std::string& encodingScheme)
            :
            utl::xcpt(mkMsg(encodingScheme))
        {}
    };

    class BindingSiteAlreadyBoundXcpt : public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& molName, const std::string& bindingSiteName)
        {
            std::ostringstream msgStream;
            msgStream << "Mol '" 
                      << molName 
                      << "' reqested already-bound binding site '"
                      << bindingSiteName 
                      << "' to be bound.  Semantic error.";
            return msgStream.str();
        }

    public:
        BindingSiteAlreadyBoundXcpt(const std::string& molName, const std::string& bindingSiteName)
            :
            utl::xcpt(mkMsg(molName, bindingSiteName))
        {}
    };


    class BindingSiteAlreadyUnboundXcpt : public utl::xcpt
    {
        static std::string
        mkMsg(const std::string& molName, const std::string& bindingSiteName)
        {
            std::ostringstream msgStream;
            msgStream << "Mol '" 
                      << molName 
                      << "' reqested already non-bound binding site '"
                      << bindingSiteName 
                      << "' to be unbound.  Semantic error.";
            return msgStream.str();
        }

    public:
        BindingSiteAlreadyUnboundXcpt(const std::string& molName, const std::string& bindingSiteName)
            :
            utl::xcpt(mkMsg(molName, bindingSiteName))
        {}
    };
    
    DEFINE_STANDARD_MSG_EXCEPTION_CLASS( MissingNameEncoderXcpt, 
                                         "Missing Name Encoder Exception: No No ptrNameAssembler set yet!!!.");
}

#endif

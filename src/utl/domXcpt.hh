#ifndef __DOMXCPT_HH
#define __DOMXCPT_HH

#include "utl/xcpt.hh"
#include "utl/dom.hh"

// This file includes all the definitions for the member functions of dom 
// exceptions:

// badChildCountXcpt
// badDoubleAttrXcpt
// badElementCastXcpt
// badIntAttrXcpt
// badNNIntAttrXcpt
// badNNDoubleAttrXcpt
// badPosIntAttrXcpt

namespace utl
{
    namespace dom
    {

        class badNNIntAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pOffendingElt,
                   const std::string& rAttrName,
                   int badAttrValue);

        public:
            badNNIntAttrXcpt (const xmlpp::Element* pOffendingElt,
                              const std::string& rAttrName,
                              int badAttrValue);
        };

        class badNNDoubleAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pOffendingElement,
                   const std::string& rAttrName,
                   double badDoubleValue);

        public:
            badNNDoubleAttrXcpt (const xmlpp::Element* pOffendingElement,
                                 const std::string& rAttrName,
                                 double badDoubleValue)
                throw();
        };


        class badChildCountXcpt :
            public xcpt
        {
            static std::string
            mkGeneralMsg (const xmlpp::Node* pParentNode,
                          const std::string& rChildName,
                          int requiredCount,
                          int actualCount);

            static std::string
            mkChoiceMsg (const xmlpp::Node* pParentNode,
                         int actualCount);

            static std::string
            mkOneOrMoreMsg (const xmlpp::Node* pParentNode);

            // In support of RNG's "optional" construct.
            static std::string
            mkZeroOrOneMsg (const xmlpp::Node* pParentNode,
                            const std::string& rChildName,
                            int actualCount);

            // This private constructor arranges for creating messages
            // in different ways, returning the same class of exception,
            // under different circumstances.
            badChildCountXcpt (const std::string& rMsg);

        public:

            // For when a definite number of children with a particular name
            // (e.g. 1) is required, but another number appears.
            static
            badChildCountXcpt
            general (const xmlpp::Node* pParentNode,
                     const std::string& rChildName,
                     int requiredCount,
                     int actualCount)
                throw();

            // When there are several possibilities for the child element's name, but
            // only one must appear, as in an RNG schema "choice" construct.
            static
            badChildCountXcpt
            choice (const xmlpp::Node* pParentNode,
                    int actualCount)
                throw();

            // When there are several possibilities for the child element's name, and
            // at least one must appear, as in an RNG schema "oneOrMore" construct.
            static
            badChildCountXcpt
            oneOrMore (const xmlpp::Node* pParentNode)
                throw();

            // Using differently named static functions avoids the problem of
            // constructors all needing different signatures.  When it's time to use
            // it, this looks like throw badChildCountXcpt::optional(...).
            static
            badChildCountXcpt
            zeroOrOne (const xmlpp::Node* pParentNode,
                       const std::string& rChildName,
                       int actualCount)
                throw();
        };

        class badDoubleAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pElement,
                   const std::string& rAttrName,
                   const std::string& rBadAttrValue);

        public:
            badDoubleAttrXcpt (const xmlpp::Element* pElement,
                               const std::string& rAttrName,
                               const std::string& rBadAttrValue)
                throw();
        };

// Thrown when a node is unexpectedly not an element.
        class badElementCastXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Node* pNode);

        public:
            badElementCastXcpt (const xmlpp::Node* pNode)
                throw();
        };

        class badIntAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pOffendingElement,
                   const std::string& rAttrName,
                   const std::string& rBadAttrValue);

        public:
            badIntAttrXcpt (const xmlpp::Element* pOffendingElement,
                            const std::string& rAttrName,
                            const std::string& rBadAttrValue);
        };

        class badPosDoubleAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pOffendingElement,
                   const std::string& rAttrName,
                   double badDoubleValue);

        public:
            badPosDoubleAttrXcpt (const xmlpp::Element* pOffendingElement,
                                  const std::string& rAttrName,
                                  double badDoubleValue)
                throw();
        };

        class badPosIntAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pOffendingElt,
                   const std::string& rAttrName,
                   int badAttrValue);

        public:
            badPosIntAttrXcpt (const xmlpp::Element* pOffendingElt,
                               const std::string& rAttrName,
                               int badAttrValue);
        };

        class missingAttrXcpt :
            public xcpt
        {
            static std::string
            mkMsg (const xmlpp::Element* pElt,
                   const std::string& rAttrName);

        public:
            missingAttrXcpt (const xmlpp::Element* pElt,
                             const std::string& rAttrName)
                throw();
        };

        // Used in domBatchJob.
        class noDocumentParsedXcpt :
                    public xcpt
        {
        public:
            noDocumentParsedXcpt (void) :
                    xcpt ("Test of parser shows no document has been parsed.")
            {}
        };
    }
}


#endif

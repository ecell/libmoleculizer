#ifndef RXNNETWORKCATALOG_HH
#define RXNNETWORKCATALOG_HH

#include "utl/autoCatalog.hh"
#include <vector>
#include <list>
#include <stack>

namespace mzr
{

    // Both speciesType and reactionType must have a public function
    // expand reaction network, ie each derived from ReactionNetworkComponent.

    template <typename speciesType,
              typename reactionType>
    class ReactionNetworkDescription
    {
    
    public:
        typedef speciesType SpeciesType;
        typedef SpeciesType* ptrSpeciesType;
    
        typedef reactionType ReactionType;
        typedef ReactionType* ptrReactionType;

        typedef utl::catalog<SpeciesType> SpeciesCatalog;
        typedef utl::catalog<ReactionType> ReactionCatalog;
        typedef std::list<ptrSpeciesType> SpeciesList;
        typedef std::list<ptrReactionType> ReactionList;
    
        typedef typename SpeciesList::iterator SpeciesListIter;
        typedef typename SpeciesList::const_iterator SpeciesListCIter;
        typedef SpeciesListCIter SpeciesToken;
    
        typedef typename ReactionList::iterator ReactionListIter;
        typedef typename ReactionList::const_iterator ReactionListCIter;
        typedef ReactionListCIter ReactionToken;
    
        typedef std::stack<SpeciesListCIter> SpeciesDeltaStack;
        typedef std::stack<ReactionListCIter> ReactionDeltaStack;

    public:

        const SpeciesCatalog& 
        getSpeciesCatalog() const
        {
            return (*theSpeciesCatalog);
        }

        const ReactionCatalog&
        getReactionCatalog() const
        {
            return *theReactionCatalog;
        }


        // Doesn't memory manage anything
        bool 
        recordSpecies( ptrSpeciesType pSpecies)
        {
            const bool speciesNewToCatalog = theSpeciesCatalog->addEntry( pSpecies->getName(),
                                                                          pSpecies);
            if (speciesNewToCatalog) 
            {
                theSpeciesList->push_back( pSpecies );
                ++totalSpeciesCreated;
                ++deltaSpeciesCreated;
            }

            return speciesNewToCatalog;
        }

        bool
        recordReaction( ptrReactionType pRxn )
        {
            // Ensure that "getCanonicalName" is correct.
            const bool rxnNewToCatalog = theReactionCatalog->addEntry( pRxn->getCanonicalName(), 
                                                                       pRxn );
            if (rxnNewToCatalog) 
            {
                theReactionList->push_back( pRxn );
                ++totalSpeciesCreated;
                ++deltaSpeciesCreated;
            }

            return rxnNewToCatalog;
        }

        void mustRecordSpecies( ptrSpeciesType pSpecies ) 
            throw( std::exception )
        {
            if ( !recordSpecies( pSpecies ) )
            {
                // Do something terrible.
            }
        }
  
        void
        mustRecordReaction( ptrReactionType pRxn) 
            throw( std::exception )
        {
            if ( !recordReactions( pRxn ) )
            {
                // Do something terrible
            }
        }

        int 
        getTotalNumberReactions() const
        {
            return totalReactionsCreated;
        }

        int 
        getTotalNumerSpecies() const
        {
            return totalSpeciesCreated;
        }

        int 
        getDeltaNumberReactions() const
        {
            return deltaReactionsCreated;
        }

        int 
        getDeltaNumberSpecies() const
        {
            return deltaSpeciesCreated;
        }

        const SpeciesDeltaStack&
        getSpeciesDeltaStack() const
        {
            return theSpeciesDeltaStack;
        }

        const ReactionDeltaStack&
        getReactionDeltaStack() const
        {
            return theReactionDeltaStack;
        }

        // Is this ok?
        const ReactionListCIter&
        getLastReactionNetworkCreationDelta() const
        {
            return theReactionDeltaStack.top();
        }

        // Is this ok?
        SpeciesListCIter 
        getLastSpeciesCreationDelta() const
        {
            return theSpeciesDeltaStack.top();
        }
  
        void recordCurrentState()
        {
            deltaSpeciesCreated = 0;
            deltaReactionsCreated = 0;

            theSpeciesDeltaStack.push( theSpeciesList.end() );
            theReactionDeltaStack.push( theSpeciesList.end() );
        }
    
        void incrementNetworkBySpeciesName(const std::string& rName)
        {
            ptrSpeciesType ptrSpecies = theSpeciesCatalog->findEntry(rName);
            if ( !ptrSpecies ) return;

            recordCurrentState();
            ptrSpecies->expandReactionNetwork();
        }

        void incrementNetworkBySpeciesToken(const SpeciesToken& specTok)
        {
            // Is there way to check this is a valid iterator to the SpeciesList?
            // Probably not....

            recordCurrentState();
            (*specTok)->expandReactionNetwork();
        }

        void incrementNetworkByReactionName(const std::string& rName)
        {
            ptrReactionType pRxn = theReactionCatalog->findEntry( rName );
      
            if (!pRxn) return;

            recordCurrentState();
            theReactionCatalog->findEntry( rName )->expandReactionNetwork();
        }
    
        void incrementNetworkByReactionToken(const SpeciesToken& rxnTok)
        {
            recordCurrentState();
            (*rxnTok)->expandReactionNetwork();
        }

        ~ReactionNetworkDescription()
        {}

        ReactionNetworkDescription()
            : 
            theSpeciesCatalog( NULL ),
            theSpeciesList( NULL ),
            theReactionCatalog( NULL ),
            theReactionList( NULL )
        {
        }

        ReactionNetworkDescription( SpeciesCatalog* aSpeciesCatalog,
                                    SpeciesList* aSpeciesList,
                                    ReactionCatalog* aReactionCatalog,
                                    ReactionList* aReactionList)
            :
            theSpeciesCatalog( aSpeciesCatalog ),
            theSpeciesList( aSpeciesList ),
            theReactionCatalog( aReactionCatalog ),
            theReactionList( aReactionList ),
            totalSpeciesCreated( 0 ),
            totalReactionsCreated( 0 )
        {
            recordCurrentState();
        }


        void configureDataRepository( SpeciesCatalog* aSpeciesCatalog,
                                      SpeciesList* aSpeciesList,
                                      ReactionCatalog* aReactionCatalog,
                                      ReactionList* aReactionList )
        {
            theSpeciesCatalog = aSpeciesCatalog;
            theSpeciesList = aSpeciesList;
            theReactionCatalog = aReactionCatalog;
            theReactionList = aReactionList;

            // I don't know exactly what I should do here,
            // but I think I should just reset everything.

        }

    private:

        SpeciesCatalog* theSpeciesCatalog;
        SpeciesList*    theSpeciesList;

        ReactionCatalog* theReactionCatalog;
        ReactionList*    theReactionList;

        // Could these get too big?...
        SpeciesDeltaStack theSpeciesDeltaStack;
        ReactionDeltaStack theReactionDeltaStack;

        unsigned int deltaSpeciesCreated;
        unsigned int deltaReactionsCreated;

        unsigned int totalSpeciesCreated;
        unsigned int totalReactionsCreated;


    };    

    
}
#endif // RXNNETWORKCATALOG_HH
    

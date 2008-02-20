#ifndef RXNNETWORKCATALOG_HH
#define RXNNETWORKCATALOG_HH

#include "utl/autoCatalog.hh"
#include "utl/xcpt.hh"
#include <vector>
#include <sstream>
#include <list>
#include <stack>

namespace mzr
{

    class DuplicatedCatalogEntryXcpt : public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& ElementName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: "
                      << "Object " 
                      << "'"
                      << ElementName 
                      << "' was already present in the ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
        DuplicatedCatalogEntryXcpt( const std::string& refObjectName) 
            :
            utl::xcpt( mkMsg( refObjectName ) )
        {}
    };

    class NoSuchSpeciesXcpt : public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& speciesName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: Species '"
                      << speciesName 
                      << "' was not found in ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
        NoSuchSpeciesXcpt( const std::string& speciesName )
            :
            xcpt( mkMsg( speciesName ) )
        {}
        
    };

    class NoSuchReactionXcpt : public utl::xcpt
    {
        static std::string
        mkMsg( const std::string& reactionName )
        {
            std::ostringstream msgStream;
            msgStream << "Internal Error: Reaction '"
                      << reactionName 
                      << "' was not found in ReactionNetworkDescription.";
            return msgStream.str();
        }

    public:
        NoSuchReactionXcpt( const std::string& reactionName )
            :
            xcpt( mkMsg( reactionName ) )
        {}
    };

            

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

        // The catalogs record everything in the reaction network.
        // There are two of them, the SpeciesCatalog and the ReactionCatalog.
        const SpeciesCatalog& 
        getSpeciesCatalog() const
        {
            return (*theSpeciesCatalog);
        }

        const ReactionCatalog&
        getReactionCatalog() const
        {
            return (*theReactionCatalog);
        }

        // These two functions are the interface that reaction generators use 
        // to record their species.  They add the (speciesName, speciesPtr)
        // and (reactionName, reactionPtr) entries to the catalog respectively.

        // If the object being recorded is new, we also record it as a hit in the
        // total number of species/reactions as well as the delta number of species/
        // reactions.
        
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
                ++totalReactionsCreated;
                ++deltaReactionsCreated;
            }

            return rxnNewToCatalog;
        }


        
        // These are not used by anyone at the moment.  It might be better to use them instead
        // of the plain old recordSpecies/recordReaction functions, but who knows.
        void mustRecordSpecies( ptrSpeciesType pSpecies ) throw(utl::xcpt)
        {
            if ( !recordSpecies( pSpecies ) )
            {
                throw DuplicatedCatalogEntryXcpt( pSpecies->getName() );
            }
        }
  
        void
        mustRecordReaction( ptrReactionType pRxn) throw(utl::xcpt)
        {
            if ( !recordReactions( pRxn ) )
            {
                throw DuplicatedCatalogEntryXcpt( pRxn->getName() );
            }
        }

        // The getTotalNumber functions are fairly obvious.  Every time a new species or reaction
        // is recorded, the member variables totalNumberSpecies or totalNumberReactions is 
        // incremented.  This number is returned by these two functions.  
        unsigned int 
        getTotalNumberSpecies() const
        {
            return totalSpeciesCreated;
        }

        unsigned int 
        getTotalNumberReactions() const
        {
            return totalReactionsCreated;
        }

        // The getDeltaNumber(Species|Reactions) functions return the number
        // of those objects created since the last time the recordCurrentState
        // function was called.
        unsigned int 
        getDeltaNumberReactions() const
        {
            return deltaReactionsCreated;
        }

        unsigned int 
        getDeltaNumberSpecies() const
        {
            return deltaSpeciesCreated;
        }


        // This returns a stack of iterators to the speciesList and reactionList.
        // Both work according to the following principle.

        // Bottom| I_0, I_1, I_2, ...., I_n | top
        // Such that the range [I_j, I_k) where j<k, will consist of all the 
        // (species|reactions) that were created between the jth and kth callings 
        // of "recordCurrentState".  Among other things 
        // [get(Species|Reaction)DeltaStack.top(), the(Species|Reactions)List.end() )
        // will list all (species|reactions) that have been recorded since the last time 
        // getCurrentState was called.

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

        // THese two are self-explanatory.  They were created in order to easily be able to 
        // iterate through [getLastSpeciesCreationDelta(), theSpeciesList.end() ) and 
        // [getLastReactionCreationDelta(), theReactionList.end() ), which will give all the 
        // (species|reactions) created since the last time recordCurrent state was called.
        const SpeciesListCIter&
        getLastSpeciesCreationDelta() const
        {
            return theSpeciesDeltaStack.top();
        }

        const ReactionListCIter&
        getLastReactionCreationDelta() const
        {
            return theReactionDeltaStack.top();
        }

        // This is QUITE IMPORTANT.  This function is used as the reset button here.  
        // This interface is only told to record new species and reactions.  This 
        // recording typically happens automatically as a consequence of species/reaction
        // generation.  Thus the user needs to call this function when appropriate.

        // [ SIDEBAR -- It could be possible to have the ReactionNetworkComponents   ]
        // [ SIDEBAR -- call this right before ending their expandNetwork member     ]
        // [ SIDEBAR -- functions.  NJA                                              ]

        // What this function actually does is reset the deltaNumberCreated indexes and 
        // pushes the back end of the NewSpecies/NewReactions lists onto a thing
        // caled the(Species|Reaction)DeltaStack.
        void recordCurrentState()
        {
            deltaSpeciesCreated = 0;
            deltaReactionsCreated = 0;

            theSpeciesDeltaStack.push( theSpeciesList->end() );
            theReactionDeltaStack.push( theReactionList->end() );
        }



        // These four functions are largely surperfluous, just nice interface 
        // functions for
        void incrementNetworkBySpeciesName(const std::string& rName) throw(utl::xcpt)
        {
            ptrSpeciesType ptrSpecies = theSpeciesCatalog->findEntry(rName);

            // Throw an exception if we can't find our man.
            if ( !ptrSpecies ) throw NoSuchSpeciesXcpt( rName );

            recordCurrentState();
            ptrSpecies->expandReactionNetwork();
        }

        void incrementNetworkByReactionName(const std::string& rName) throw(utl::xcpt)
        {
            ptrReactionType pRxn = theReactionCatalog->findEntry( rName );
      
            // Throw an exception if we can't find our man.
            if (!pRxn) throw NoSuchReactionXcpt( rName );

            recordCurrentState();
            pRxn->expandReactionNetwork();
        }


//         void incrementNetworkBySpeciesToken(const SpeciesToken& specTok)
//         {
//             // Is there way to check this is a valid iterator to the SpeciesList?
//             // Probably not....

//             recordCurrentState();
//             (*specTok)->expandReactionNetwork();
//         }


    
//         void incrementNetworkByReactionToken(const SpeciesToken& rxnTok)
//         {
//             recordCurrentState();
//             (*rxnTok)->expandReactionNetwork();
//         }

        ~ReactionNetworkDescription()
        {

            // DO NOT DELETE THESE OBJECTS.  Because they are provided by an inheriting class, 
            // they probably don't even exist at the point of this function.

            // DO NOT delete theSpeciesCatalog;
            // DO NOT delete theSpeciesList;
            // DO NOT delete theReactionCatalog;
            // DO NOT delete theReactionList;

            // But we DO memory manage SpeciesDeltaStack and ReactionDeltaStack.
            // This is currently allocated on the stack, and so is automatically 
            // memory managed.  I suspect this may change when memory usage is profiled.
        }

        ReactionNetworkDescription()
            : 
            theSpeciesCatalog( NULL ),
            theSpeciesList( NULL ),
            theReactionCatalog( NULL ),
            theReactionList( NULL ),
            totalSpeciesCreated( 0 ),
            totalReactionsCreated( 0 )
        {}

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


        // This is kinda bogus and I don't much care for it.  The problem is that
        // This class is made to be inherited from. (in this application this is done by the 
        // moleculizer class).  Unfortunately I would like the inheriting class to provide its 
        // data storage here (although this class does manage the DeltaStacks of iterators), and
        // so this class will usually be constructed prior to the storage provided by its child.
        // As near as I can tell, the implication is that this may be the best way I can go.  
        void configureDataRepository( SpeciesCatalog* aSpeciesCatalog,
                                      SpeciesList* aSpeciesList,
                                      ReactionCatalog* aReactionCatalog,
                                      ReactionList* aReactionList )
        {
            theSpeciesCatalog = aSpeciesCatalog;
            theSpeciesList = aSpeciesList;
            theReactionCatalog = aReactionCatalog;
            theReactionList = aReactionList;

            
            // We decide to be agnostic about the data located here; however, we must blow
            // out the deltaStacks because we can't risk having bad iterators to non-existant
            // data.  
            clearDeltaStacks();
            recordCurrentState();
        }

    private:

        void clearDeltaStacks()
        {
            this->clearStack( theSpeciesDeltaStack );
            this->clearStack( theReactionDeltaStack );
        }

        
        template <typename T>
        void clearStack( std::stack<T>& aStack)
        {
            while (! aStack.empty() )
            {
                aStack.pop();
            }
        }

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
    

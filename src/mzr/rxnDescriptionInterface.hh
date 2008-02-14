#ifndef RXNDESCRIPTIONINTERFACE_HH
#define RXNDESCRIPTIONINTERFACE_HHo

#include "utl/autoCatalog.hh"
#include <vector>
#include <list>

namespace mzr
{
    class mzrReaction;
    class mzsSpecies;
}

struct ReactionNetworkDescription
{
  // Doesn't memory manage anything

  bool 
  recordSpecies( mzrSpecies* pSpecies)
  {}

  bool
  recordReaction( mzrReaction* pReaction)
  {
    std::string reactionName = pReaction->getCanoncalName();
    return reactionCatalog->addEntry( reactionName, pReaction );
  }

  utl::catalog<mzr::mzrSpecies>* pSpeciesCatalog;
  std::list<mzr::mzrSpecies*>* pSpeciesList;
  utl::catalog<mzr::mzrReaction>* pReactionCatalog;
  std::list<mzr::mzrReaction*>* pReactionList;

  ReactionNetowrkDescription( utl::catalog<mzr::mzrSpecies>* aSpeciesCatalog,
			      std::list<mzr::mzrSpecies*>* aSpeciesList,
			      utl::catalog<mzr::mzrReaction>* aReactionCatalog,
			      std::list<mzr::mzrReaction*>* aReactionList) 
    :
    pSpeciesCatalog( aSpeciesCatalog ),
    pSpeciesList( aSpeciesList ),
    pReactionCatalog( aReactionCatalog ),
    pReactionList( aReactionList )
  {}

};


#endif

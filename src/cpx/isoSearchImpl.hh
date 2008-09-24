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


#ifndef CPX_ISOSEARCHIMPL_H
#define CPX_ISOSEARCHIMPL_H

#include "cpx/plexNotConnectedXcpt.hh"
#include "nauty/nauty.h"

namespace cpx
{
template<class plexT>
bool
isoSearch<plexT>::
findIso(void) const
{
if((rLeft.bindings.size() != rRight.bindings.size())
||
(rLeft.mols.size() != rRight.mols.size()))
{
return false;
}
else
{
return nautyIsomorphismCheck(rLeft, rRight);
}

}


template<class plexT>
void
isoSearch<plexT>::enforcePlexIsSimpleGraph( const plexT& plex) const
throw(plexIsNotSimpleGraphXcpt)
{
typedef std::pair<molSpec, molSpec> KeyType;
typedef int ValueType;

std::map< KeyType, ValueType> edgeMap;

BOOST_FOREACH(const binding& bnd, plex.bindings)
{
molSpec leftMolSpec = bnd.leftSite().molNdx();
molSpec rightMolSpec = bnd.rightSite().molNdx();

std::pair<molSpec, molSpec> molPair;
rightMolSpec < leftMolSpec? molPair = std::make_pair( rightMolSpec, leftMolSpec) :
molPair = std::make_pair(leftMolSpec, rightMolSpec);

// Checks to see if the binding is intra-molecular.
if (leftMolSpec == rightMolSpec)
throw plexIsNotSimpleGraphXcpt();

std::pair< std::map<KeyType, ValueType>::iterator, bool> insertResult
= edgeMap.insert( std::make_pair( molPair, 1));

if (!insertResult.second) throw plexIsNotSimpleGraphXcpt();
}

return;
}


template <class plexT>
bool
isoSearch<plexT>::nautyIsomorphismCheck(const plexT& rLeftPlex,
const plexT& rRightPlex) const
{
// CREATING OPTIONS
// This should ultimately be moved elsewhere, as in a static member of the class.
DEFAULTOPTIONS_GRAPH(options);
options.getcanon = TRUE;
options.defaultptn = FALSE;
statsblk stats;


// ENSURING THAT THE SAME NODE SIZE WILL WORK FOR BOTH COMPLEXES
// Semantically this is unnecessary, as this check has/should have already been performed
// by any code that calls this function. However, within this function it guarantees that n
// is the same between both complexes.
if (rLeftPlex.mols.size() != rRightPlex.mols.size()) return false;
const int complexSize( rLeftPlex.mols.size() );


int m = (complexSize + WORDSIZE - 1) / WORDSIZE;

char error[] = "malloc";

DYNALLSTAT(graph,graph1,g1_sz);
DYNALLOC2(graph,graph1,g1_sz,m,complexSize,error);

DYNALLSTAT(graph,graph2,g2_sz);
DYNALLOC2(graph,graph2,g2_sz,m,complexSize,error);

DYNALLSTAT(graph, canong1, canong1_sz);
DYNALLSTAT(graph, canong2, canong2_sz);
DYNALLOC2(graph, canong1, canong1_sz, m, complexSize, error);
DYNALLOC2(graph, canong2, canong2_sz, m, complexSize, error);

DYNALLSTAT(int,orbits1,orbits1_sz);
DYNALLOC1(int, orbits1, orbits1_sz, complexSize, error);

DYNALLSTAT(int,orbits2,orbits2_sz);
DYNALLOC1(int, orbits2, orbits2_sz, complexSize, error);

DYNALLSTAT(setword,workspace,workspace_sz);
DYNALLOC1(setword, workspace, workspace_sz, 100*m, error);

std::vector<int> labelingMap1;
std::vector<int> coloringPartition1;

createGraphFromPlex(rLeftPlex, graph1, labelingMap1, coloringPartition1, m );

std::vector<int> labelingMap2;
std::vector<int> coloringPartition2;

createGraphFromPlex(rRightPlex, graph2, labelingMap2, coloringPartition2, m );

// COMPARE TO MAKE SURE THAT THE COLORINGS ARE THE SAME HERE.....
if (!compareColoringPartitions(rLeftPlex, labelingMap1, coloringPartition1,
rRightPlex, labelingMap2, coloringPartition2) )
{
return false;
}

// CALL NAUTY ON EVERYTHING AND COMPARE THEIR CANONICAL GRAPHS TO DETERMINE WHETHER THEY ARE ISOMORPHIC.
nauty(graph1, &labelingMap1[0], &coloringPartition1[0], NULL, orbits1, &options, &stats, workspace, 100*m, m, complexSize, canong1);
nauty(graph2, &labelingMap2[0], &coloringPartition2[0], NULL, orbits2, &options, &stats, workspace, 100*m, m, complexSize, canong2);
const bool areIsomorphic( *canong1 == *canong2);

// IF ISOMORPHIC CONSTRUCT A MAPPING BETWEEN THEM AND CALL 'onSuccess'.
if (areIsomorphic)
{
// Use the maps to create a ple
plexIso tmpIso(rLeftPlex.mols.size(),
rLeftPlex.bindings.size(),
rRightPlex.mols.size(),
rRightPlex.bindings.size());

createMappingBetweenIsomorphicPlexes(tmpIso,
rLeftPlex,  labelingMap1,
rRightPlex, labelingMap2);

onSuccess(tmpIso);
}

// DELETE ALL MEMORY ALLOCATED BY THIS FUNCTION.
// todo.

return areIsomorphic;

}

template<class plexT>
void isoSearch<plexT>::createGraphFromPlex( const plexT& aPlex,
graph* theGraph,
std::vector<int>& labelingMap,
std::vector<int>& ptnMap,
int m) const
{
// The purpose of this function is to take aPlex and use it to initialize the
// fundamental nauty data structures, a graph, a labeling, and a partitionMap.
// NOTE -- I assume that theGraph is properly initialized to the correct size
//         outside of the function.  I don't know if this is the best way to do
//         things, but if this condition is not met, things will go downhill fast.
//
// Graph -- this is the particular form of graph that Nauty likes, with its own
//          idiosyncratic, C + Macro way of constructing it.
//
// Labeling -- This is a permutation on the edge nodes of the graph such that the
//             resulting node list has edges of the same color class adjacent.
//
// Partition -- This is an array whose values are 1 everywhere except where it's 0. :-)
//              A value is 0 at an index 'ndx' if edge[labeling[ndx]] is of one color
//              and edge[labeling[ndx + 1]] is of another.

const int complexSize( aPlex.mols.size() );

labelingMap.clear();
ptnMap.clear();
labelingMap.resize( complexSize, 0 );
ptnMap.resize( complexSize, 1 );


// INITIALIZING DATA STRUCTURES
std::vector<std::set<int>* > edgeMap;
for(int ii = 0; ii != (int) aPlex.mols.size(); ++ii)
{
edgeMap.push_back(new std::set<int>() );

// I am here creating the initial labeling map even though it is totally distinct
// from the edgeMap.
labelingMap[ii] = ii;
}


// CREATING THE GRAPH

BOOST_FOREACH( const binding& bnd, aPlex.bindings)
{
int molNdx1 = bnd.leftSite().molNdx();
int molNdx2 = bnd.rightSite().molNdx();

edgeMap[molNdx1]->insert(molNdx2);
edgeMap[molNdx2]->insert(molNdx1);
}


// Now to put all this information into the graph.
set* gv;
for(unsigned int nodeNdx = 0; nodeNdx != edgeMap.size(); ++nodeNdx)
{
gv = GRAPHROW(theGraph, nodeNdx, m);
EMPTYSET(gv, m);

for(std::set<int>::const_iterator iter = edgeMap[nodeNdx]->begin();
iter != edgeMap[nodeNdx]->end();
++iter)
{
ADDELEMENT(gv, *iter);
}
}

// CREATING THE LABELING MAP

// Let us create the labeling map, by using this function object which stores aPlex
// and returns ndx1 < ndx2 iff aPlex.mols[ndx1].getName() < aPlex.mols[ndx2].getName().
std::sort(labelingMap.begin(),
labelingMap.end(),
IndexArraySorter(aPlex));


// CREATING THE PARTITION MAP
// First off, the last element will always be 0, because it always ends a color partition.
ptnMap[ ptnMap.size() - 1 ] = 0;

for(unsigned int ii = 0; ii != ptnMap.size() -1; ++ii)
{
if( aPlex.mols[ labelingMap[ii] ]->getName() != aPlex.mols[labelingMap[ii + 1] ] ->getName())
{
ptnMap[ ii ] = 0;
}
else
{
ptnMap[ii] = 1;
}
}

}


template<class plexT>
bool
isoSearch<plexT>::
mapRestBindings(int leftBindingNdx,
const plexIso& rCurrentIso) const
{
// Are we done?
if(((int) rLeft.bindings.size()) <= leftBindingNdx)
{
onSuccess(rCurrentIso);
return true;
}

// Try to extend the given isomorphism over the given left binding
// until one is found that can be extended to a full isomorpism.
for(int rightBindingNdx = 0;
rightBindingNdx < (int) rRight.bindings.size();
rightBindingNdx++)
{
// Construct a new temporary isomorphism for this trial.
plexIso tmpIso(rCurrentIso);

if(tmpIso.tryMapBinding(rLeft,
leftBindingNdx,
rRight,
rightBindingNdx)
&& mapRestBindings(leftBindingNdx + 1,
tmpIso))
{
return true;
}
}
return false;
}


template<class plexT>
bool
isoSearch<plexT>::
findInjection(void) const
{
plexIso tmpIso(rLeft.mols.size(),
rLeft.bindings.size(),
rRight.mols.size(),
rRight.bindings.size());

// Since we're assuming that the left plex is connected,
// we've got two cases, either every mol is on a binding
// or there is only one mol.
if(rLeft.bindings.size() > 0)
{
// Search for mappings of all the bindings in the pattern.
if(mapRestBindings(0,
tmpIso))
{
return true;
}
}
else if(1 == rLeft.mols.size())
{
// Attempt to match with each of the mols in the target complex.
for(int tgtMolNdx = 0;
tgtMolNdx < (int) rRight.mols.size();
tgtMolNdx++)
{
if(tmpIso.forward.canMapMol(rLeft,
0,
rRight,
tgtMolNdx))
{
tmpIso.forward.molMap[0] = tgtMolNdx;
tmpIso.backward.molMap[tgtMolNdx] = 0;
onSuccess(tmpIso);
return true;
}
}
}
else
{
throw plexNotConnectedXcpt();
}

return false;
}


template <class plexT>
bool
isoSearch<plexT>::compareColoringPartitions(const plexT& plex1, const std::vector<int>& labelingMap1, const std::vector<int>& partitionMap1,
const plexT& plex2, const std::vector<int>& labelingMap2, const std::vector<int>& partitionMap2) const
{
// Replace with throw BadGraphSpecificationXcpt("");
if (labelingMap1.size() != partitionMap1.size() ||
labelingMap2.size() != partitionMap2.size() ) throw 666;


// Couple definitions:
// molMap - a vector of the mols in the plex.
// labelMap - a permutation which sorts mols by color
// partition Map -- an array which is 1 on the last entity in the coloring partition and 0 everywhere else

// Example molMap -- [D, A, B, C, B]
// labeling map   -- [4, 0, 1, 3, 2]
// Partition map  -- [1, 0, 1, 1, 1]

if( partitionMap1 != partitionMap2 ) return false;

const std::vector<int>& partitionMap( partitionMap1 );


std::vector<int>::const_iterator begin_first = labelingMap1.begin();
std::vector<int>::const_iterator begin_second = labelingMap2.begin();

for(unsigned int partNdx = 0; partNdx != partitionMap.size(); ++partNdx)
{
if( partitionMap[partNdx] == 0) continue;

std::vector<int>::const_iterator indexIter_1 = find(labelingMap1.begin(), labelingMap1.end(), partNdx);
if (indexIter_1 == labelingMap1.end())
{
std::cerr << "Error in isoSearchImpl.hh(1)" << endl;
throw 666;
}

std::vector<int>::const_iterator indexIter_2 = find(labelingMap2.begin(), labelingMap2.end(), partNdx);

if (indexIter_2 == labelingMap2.end())
{
std::cerr << "Error in isoSearchImpl.hh(2)" << endl;
throw 666;
}

int index_1 = indexIter_1 - begin_first;
int index_2 = indexIter_2 - begin_second;

if (plex1.mols[index_1]->getName() != plex2.mols[index_2]->getName()) return false;

}

return true;
}

template <class plexT>
void
isoSearch<plexT>::createMappingBetweenIsomorphicPlexes(plexIso& isoMapping,
const plexT& leftPlex,  std::vector<int>& leftCanonicalLabeling,
const plexT& rightPlex, std::vector<int>& rightCanonicalLabeling) const
{

std::vector<int> leftPermutation;
std::vector<int> rightPermutation;

createPermutationFromRelabeling(leftCanonicalLabeling, leftPermutation);
createPermutationFromRelabeling(rightCanonicalLabeling, rightPermutation);

static int xxx = 0;


std::vector<int> isoToRightMapping(rightPermutation.size(), -1);
invertMapping( rightPermutation, isoToRightMapping);

// CREATE THE MOL SUBMAPPINGS
std::vector<int> leftMolsToRightMolsMapping( leftPlex.mols.size(), -1);
for(int leftMolNdx = 0; leftMolNdx != (int) leftPlex.mols.size(); ++leftMolNdx)
{
int rightMolNdx = isoToRightMapping[ leftPermutation[ leftMolNdx ]];

isoMapping.forward.molMap[leftMolNdx] = rightMolNdx;
isoMapping.backward.molMap[rightMolNdx] = leftMolNdx;

if (leftPlex.mols[leftMolNdx]->getName() != rightPlex.mols[rightMolNdx]->getName())
{
std::cerr << xxx << ":\tERR" << std::endl;
}
}

xxx++;


// CREATE THE MAPPINGS BETWEEN BINDINGS.
// This is probably slower than it needs to be, although it's merely a slow linear search.

std::vector<binding>::const_iterator rightBindingBeginIter = rightPlex.bindings.begin();

for(int leftPlexBindingNdx = 0; leftPlexBindingNdx != (int) leftPlex.bindings.size(); ++leftPlexBindingNdx)
{
// This is also a lot more bit-twiddly than it ought to be.
const binding& theBinding = leftPlex.bindings[leftPlexBindingNdx];

const siteSpec& leftSite = theBinding.leftSite();
const siteSpec& rightSite = theBinding.rightSite();

int mappedLeftMolNdx = isoMapping.forward.molMap[leftSite.molNdx()];
int mappedRightMolNdx = isoMapping.forward.molMap[rightSite.molNdx()];

binding mappedBinding( siteSpec(mappedLeftMolNdx,  leftSite.siteNdx()),
siteSpec(mappedRightMolNdx, rightSite.siteNdx()));


std::vector<binding>::const_iterator iterToBindingInRightPlex = find_if( rightPlex.bindings.begin(),
rightPlex.bindings.begin(),
bindingIsEqual(mappedBinding));

if (iterToBindingInRightPlex == rightPlex.bindings.end())
{
std::cerr << "Error finding binding that SHOULD be there in isoSearchImpl.hh" << std::endl;
throw 666;
}

int rightPlexBindingNdx = iterToBindingInRightPlex - rightBindingBeginIter;

// Add (leftPlexBindingNdx -> rightPlexBindingNdx) to the forwardBindingMap.

isoMapping.forward.bindingMap[leftPlexBindingNdx] = rightPlexBindingNdx;
isoMapping.forward.bindingMap[rightPlexBindingNdx] = leftPlexBindingNdx;
}
}


template <class plexT>
void
isoSearch<plexT>::invertMapping( const std::vector<int>& mapToInvert, std::vector<int>& invertedMapping) const
{
invertedMapping.clear();
invertedMapping.resize( mapToInvert.size() );

for(int domainElement = 0; domainElement != (int) mapToInvert.size(); ++domainElement)
{
int rangeElement = mapToInvert[domainElement];
invertedMapping[rangeElement] = domainElement;
}

return;
}


template <class plexT>
void
isoSearch<plexT>::createPermutationFromRelabeling(const std::vector<int>& labelingMap,
std::vector<int>& permutationMap) const
{
permutationMap.clear();
permutationMap.resize( labelingMap.size(), -1);

std::vector<int>::const_iterator beginIter = labelingMap.begin();
std::vector<int>::const_iterator findIter = labelingMap.end();
std::vector<int>::const_iterator endIter = labelingMap.end();

for(int domainElement = 0; domainElement != (int) labelingMap.size(); ++domainElement)
{
findIter = std::find(labelingMap.begin(),
labelingMap.end(),
domainElement);

if (findIter == endIter) throw 666;

int rangeElement = findIter - beginIter;

permutationMap[domainElement] = rangeElement;
}

if (std::count(permutationMap.begin(),
permutationMap.end(),
-1) > 0) throw 666;

return;
}



}


#endif // CPX_ISOSEARCHIMPL_H

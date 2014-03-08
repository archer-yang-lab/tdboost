//------------------------------------------------------------------------------
//  by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_factory.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   manager for allocation and destruction of all nodes 
//        	  
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODEFACTORY_H
#define NODEFACTORY_H

#include <stack>
#include <list>
#include "node_terminal.h"
#include "node_continuous.h"
#include "node_categorical.h"

#define NODEFACTORY_NODNPtweedie_RESERVE ((unsigned long)50)

using namespace std;

class CNodeFactory
{
public:
    CNodeFactory();
    ~CNodeFactory();

    NPtweedieRESULT Initialize(unsigned long cDepth);
    CNodeTerminal* GetNewNodeTerminal();
    CNodeContinuous* GetNewNodeContinuous();
    CNodeCategorical* GetNewNodeCategorical();
    NPtweedieRESULT RecycleNode(CNodeTerminal *pNode);
    NPtweedieRESULT RecycleNode(CNodeContinuous *pNode);
    NPtweedieRESULT RecycleNode(CNodeCategorical *pNode);

private:
    stack<PCNodeTerminal> TerminalStack;
    stack<PCNodeContinuous> ContinuousStack;
    stack<PCNodeCategorical> CategoricalStack;

    CNodeTerminal* pNodeTerminalTemp;
    CNodeContinuous* pNodeContinuousTemp;
    CNodeCategorical* pNodeCategoricalTemp;

    CNodeTerminal aBlockTerminal[NODEFACTORY_NODNPtweedie_RESERVE];
    CNodeContinuous aBlockContinuous[NODEFACTORY_NODNPtweedie_RESERVE];
    CNodeCategorical aBlockCategorical[NODEFACTORY_NODNPtweedie_RESERVE];
};

#endif // NODEFACTORY_H




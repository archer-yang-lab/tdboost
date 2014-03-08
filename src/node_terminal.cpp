// The code is a modified version of gbm library originally written by Greg Ridgeway. See
// 
// Ridgeway, G. (2007). Generalized boosted models: A guide to the gbm package. R pack-
// age vignette. http://cran.r-project.org/web/packages/gbm.
//------------------------------------------------------------------------------
//  by Greg Ridgeway  Copyright (C) 2003
#include "node_terminal.h"
#include "node_factory.h"

CNodeTerminal::CNodeTerminal()
{
    isTerminal = true;
}


CNodeTerminal::~CNodeTerminal()
{
    #ifdef NOISY_DEBUG
    Rprintf("terminal destructor\n");
    #endif
}

NPtweedieRESULT CNodeTerminal::Adjust
(
    unsigned long cMinObsInNode
)
{
    return NPtweedie_OK;
}

NPtweedieRESULT CNodeTerminal::ApplyShrinkage
(
    double dLambda
)
{
    NPtweedieRESULT hr = NPtweedie_OK;
    dPrediction *= dLambda;

    return hr;
}



NPtweedieRESULT CNodeTerminal::Predict
(
    CDataset *pData, 
    unsigned long iRow, 
    double &dFadj
)
{
    dFadj = dPrediction;

    return NPtweedie_OK;
}




NPtweedieRESULT CNodeTerminal::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow,
    double &dFadj
)
{
    dFadj = dPrediction;

    return NPtweedie_OK;
}



NPtweedieRESULT CNodeTerminal::PrintSubtree
(
    unsigned long cIndent
)
{
    unsigned long i = 0;
    
    for(i=0; i< cIndent; i++) Rprintf("  ");
    Rprintf("N=%f, Prediction=%f *\n",
           dTrainW,
           dPrediction);

    return NPtweedie_OK;
}


NPtweedieRESULT CNodeTerminal::GetVarRelativeInfluence
(
    double *adRelInf
)
{
    return NPtweedie_OK;
}


NPtweedieRESULT CNodeTerminal::RecycleSelf
(
    CNodeFactory *pNodeFactory
)
{
    pNodeFactory->RecycleNode(this);
    return NPtweedie_OK;
};



NPtweedieRESULT CNodeTerminal::TransferTreeToRList
(
    int &iNodeID,
    CDataset *pData,
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld,
    double dShrinkage
)
{
    NPtweedieRESULT hr = NPtweedie_OK;

    aiSplitVar[iNodeID] = -1;
    adSplitPoint[iNodeID] = dShrinkage*dPrediction;
    aiLeftNode[iNodeID] = -1;
    aiRightNode[iNodeID] = -1;
    aiMissingNode[iNodeID] = -1;
    adErrorReduction[iNodeID] = 0.0;
    adWeight[iNodeID] = dTrainW;
    adPred[iNodeID] = dShrinkage*dPrediction;

    iNodeID++;

    return hr;
}



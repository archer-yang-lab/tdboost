//  NPtweedie by Yi Yang and Hui Zou  Copyright (C) 2012
#include "EDM.h"

CEDM::CEDM(double dAlpha)
{
    this->dAlpha = dAlpha;
}

CEDM::~CEDM()
{

}

// compute gradient function
NPtweedieRESULT CEDM::ComputeWorkingResponse 
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
    unsigned long nTrain
)
{
	NPtweedieRESULT hr = NPtweedie_OK;
    unsigned long i = 0;
    double dF = 0.0;
    
    if((adY == NULL) || (adF == NULL) || (adZ == NULL) || (adWeight == NULL))
    {
        hr = NPtweedie_INVALIDARG;
        goto Error;
    }
    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
			adZ[i] = -adY[i] * exp((1.0-dAlpha)*adF[i]) + exp((2.0-dAlpha)*adF[i]);
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
			dF = adF[i] + adOffset[i];
			adZ[i] = -adY[i] * exp((1.0-dAlpha)*dF) + exp((2.0-dAlpha)*dF);
        }
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




// compute likelihood function
double CEDM::Deviance
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double *adF,
    unsigned long cLength
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;
    double dF = 0.0;
    
	if(adOffset == NULL)
    {
      	for(i=0; i<cLength; i++)
      	{
         	dL += adWeight[i] * (-adY[i] * exp((1.0-dAlpha) * adF[i]) / (1.0-dAlpha) + exp((2.0-dAlpha)* adF[i]) / (2.0-dAlpha));
         	dW += adWeight[i];
      	}
    }
	else
	{
      	for(i=0; i<cLength; i++)
      	{
         	dF = adF[i] + adOffset[i];
         	dL += adWeight[i] * (-adY[i] * exp((1.0-dAlpha) * dF) / (1.0-dAlpha) + exp((2.0-dAlpha)* dF) / (2.0-dAlpha));
         	dW += adWeight[i];
      	}
    }
	
    return dL/dW;
}





// compute Initial value
NPtweedieRESULT CEDM::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
	    unsigned long i=0;
   	    double dTemp = 0.0;

        // Newton method for solving for F
        // should take about 3-6 iterations.
        double dNum=0.0;         // numerator
        double dDen=0.0;         // denominator
        double dNewtonStep=1.0;  // change
        dInitF = 0.0;
//        while(fabs(dNewtonStep) > 0.0001)
//        {
            dNum=0.0;
            dDen=0.0;
            for(i=0; i<cLength; i++)
            {
				dTemp = dInitF + ((adOffset==NULL) ? 0.0 : adOffset[i]);
                dNum += adWeight[i]*(-adY[i] * exp((1.0-dAlpha)*dTemp) + exp((2.0-dAlpha)*dTemp));
                dDen += adWeight[i]*(-adY[i] * (1.0-dAlpha) * exp((1.0-dAlpha)*dTemp) + (2.0-dAlpha) * exp((2.0-dAlpha)*dTemp));
            }
            dNewtonStep = -dNum/dDen;
            dInitF += dNewtonStep;
//        }

	//	printf("%f\n", dInitF);
        return NPtweedie_OK;

}






NPtweedieRESULT CEDM::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj
)
{

    NPtweedieRESULT hr = NPtweedie_OK;
    
    double dF = 0.0;
    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vector<double> vecdNum;	
    vector<double> vecdDen;	
    vecdNum.resize(cTermNodes);
    vecdNum.assign(vecdNum.size(),0.0);
    vecdDen.resize(cTermNodes);
    vecdDen.assign(vecdDen.size(),0.0);

	for(iObs=0; iObs<nTrain; iObs++)
	{
		if(afInBag[iObs])
		{
			dF = adF[iObs] + ((adOffset==NULL) ? 0.0 : adOffset[iObs]);
           	vecdNum[aiNodeAssign[iObs]] += 
				adW[iObs]*(-adY[iObs] * exp((1.0-dAlpha)*dF) + 
				exp((2.0-dAlpha)*dF));
           	vecdDen[aiNodeAssign[iObs]] += 
				adW[iObs]*(-adY[iObs] * (1.0-dAlpha) * exp((1.0-dAlpha)*dF) + 
				(2.0-dAlpha) * exp((2.0-dAlpha)*dF));
       	}
	}

	for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]!=NULL)
        {
			if(vecdDen[iNode] == 0)
            {
                vecpTermNodes[iNode]->dPrediction = 0.0;
            }
            else
            {
                vecpTermNodes[iNode]->dPrediction -= vecdNum[iNode]/vecdDen[iNode];
			}
        }
    }

    return hr;
}

// compute likelihood improvement after updates
double CEDM::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    double dL = 0.0;
    double dLadj = 0.0;
	double dW = 0.0;
    unsigned long i = 0;
    double dF = 0.0;
	double ddF = 0.0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
			dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
			ddF = dF + dStepSize*adFadj[i];
         	dL += adWeight[i] * (-adY[i] * exp((1.0-dAlpha) * dF) / (1.0-dAlpha) + exp((2.0-dAlpha)* dF) / (2.0-dAlpha));
         	dLadj += adWeight[i] * (-adY[i] * exp((1.0-dAlpha) * ddF) / (1.0-dAlpha) + exp((2.0-dAlpha)* ddF) / (2.0-dAlpha));    
            dW += adWeight[i];
        }
    }

    return (dL-dLadj)/dW;
}

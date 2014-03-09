//  NPtweedie by Yi Yang and Hui Zou  Copyright (C) 2012
#include "dispexpo.h"

CDispexpo::CDispexpo(double dAlpha)
{
    this->dAlpha = dAlpha;
}

CDispexpo::~CDispexpo()
{

}

// compute gradient function
NPtweedieRESULT CDispexpo::ComputeWorkingResponse 
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
double CDispexpo::Deviance
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
         	dL += adWeight[i] * (adY[i] * exp((1.0-dAlpha) * adF[i]) / (dAlpha-1.0) + exp((2.0-dAlpha)* adF[i]) / (2.0-dAlpha));
         	dW += adWeight[i];
      	}
    }
	else
	{
      	for(i=0; i<cLength; i++)
      	{
         	dF = adF[i] + adOffset[i];
         	dL += adWeight[i] * (adY[i] * exp((1.0-dAlpha) * dF) / (dAlpha-1.0) + exp((2.0-dAlpha)* dF) / (2.0-dAlpha));
         	dW += adWeight[i];
      	}
    }
	
    return dL/dW;
}






NPtweedieRESULT CDispexpo::InitF
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
                dNum += adWeight[i]*(-adY[i] * exp((1.0-dAlpha) * dTemp) + exp((2.0-dAlpha)*dTemp));
                dDen += adWeight[i]*(-adY[i] * (1.0-dAlpha) * exp((1.0-dAlpha) * dTemp) + (2.0-dAlpha) * exp((2.0-dAlpha)*dTemp));
            }
            dNewtonStep = -dNum/dDen;
            dInitF += dNewtonStep;
//        }
   
        return NPtweedie_OK;

}






NPtweedieRESULT CDispexpo::FitBestConstant
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
    
    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vector<vector<double> > vecadDiff;	
    vecadDiff.resize(cTermNodes);
	double dOffset;

    for(iObs=0; iObs<nTrain; iObs++)
    {
            if(afInBag[iObs])
            {
					dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
                    vecadDiff[aiNodeAssign[iObs]].push_back(adY[iObs]-dOffset-adF[iObs]) ;            
            }
    }
    
    
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecadDiff[iNode].size()!=0)
        {
                sort(vecadDiff[iNode].begin(),vecadDiff[iNode].end());
                double Deviance;
                double test;
                unsigned long Len = vecadDiff[iNode].size();
                unsigned long begin=0;
                unsigned long end=Len-1;

                     while (end>begin+1) 
                     {
                             Deviance=0;
                             test=floor(double(begin+end)/2.0);
                             for (unsigned long j=0; j<test; j++)
                             {
                                     Deviance += (1.0-dAlpha)*(vecadDiff[iNode][j]-vecadDiff[iNode][test]);
                             }

                             for (unsigned long j=Len-1;j>test;j--) 
                             {
                                     Deviance += dAlpha*(vecadDiff[iNode][j]-vecadDiff[iNode][test]);
                             }

                             if(Deviance>0) 
                             {
                                     begin=test;
                             }

                             else
                             {
                                     end=test;
                             }
                     }

                     double Pnum=0;
                     for (unsigned long j=0;j<(begin+1);j++) 
                     {
                             Pnum += (1.0-dAlpha)*(vecadDiff[iNode][j]);
                     }

                     for (unsigned long j=Len-1;j>(end-1);j--) 
                     {
                             Pnum += dAlpha*(vecadDiff[iNode][j]);
                     }

                     if (((1.0-dAlpha)*end + dAlpha*(Len-end))==0)
                     {
                             Pnum=0.0;  
                     }
                     else
                     {
                             Pnum=Pnum / ((1.0-dAlpha)*end + dAlpha*(Len-end));
                     }       	

                vecpTermNodes[iNode] ->dPrediction = Pnum;
                        
         }
    }	

    return hr;
}


double CDispexpo::BagImprovement
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
         	dL += adWeight[i] * (adY[i] * exp((1.0-dAlpha) * dF) / (dAlpha-1.0) + exp((2.0-dAlpha)* dF) / (2.0-dAlpha));
         	dLadj += adWeight[i] * (adY[i] * exp((1.0-dAlpha) * ddF) / (dAlpha-1.0) + exp((2.0-dAlpha)* ddF) / (2.0-dAlpha));    
            dW += adWeight[i];
        }
    }

    return (dL-dLadj)/dW;
}

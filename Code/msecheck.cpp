#include "krig.h"
double msecheck(model S)
{
        int numpts;
        int dim,dim_th;
	double mse = 0;
        MatD x0,L;
        VecD y0,theta;
	VecI th_index;
        if(!model2param_wchol(S,numpts,dim,dim_th,x0,y0,th_index,theta,L))
		std::cout<<" CHOLESKY FAILED. BACK TO WORK DUMBASS \n ";
	else
	{
//	CALCULATE MSE
		VecD y(numpts);
		predictor(S,numpts,x0,y);
		for(int i=0;i<numpts;i++)
		{
			 mse += (y0[i]-y[i])*(y0[i]-y[i]);
		}
	}
	return mse;
}

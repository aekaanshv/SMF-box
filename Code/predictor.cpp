#include "krig.h"

void predictor(model S, int num_pts1, MatD x,VecD &y)
{
// MODEL TO PARAMETERS 
	int num_pts0;
 	int dim,dim_th;
	MatD x0,L;
	VecI th_index;
 	VecD y0,theta;
	model2param(S,num_pts0,dim,dim_th,x0,y0,th_index,theta,L);


        double mean = mean_cal(num_pts0,L,y0);
	
// DUMMY VARIABLES
	VecD dummy_vec(num_pts0);
	VecD y0_minus_mean(num_pts0,0);
	int i,j;
	double sum;

// OBSERVED RESULTS - CALCULATED MEAN

	for(i=0;i<num_pts0;i++)
		y0_minus_mean[i] = y0[i] - mean;



// MULTIPLY BY INV(COVARIANCE MATRIX)

	chol_sol(num_pts0,L,y0_minus_mean,dummy_vec);

// PREDICT y : CALCULATE y = mu + (r')*(R_inv*(y0 - mu*1))

	for(i=0;i<num_pts1;i++)
	{
		sum=mean;
		VecD covx;
		expcov(num_pts0,dim,dim_th,x0,x[i],th_index,theta,covx);
		for(j=0;j<num_pts0;j++)
			sum = sum+covx[j]*dummy_vec[j];
		y[i] = sum;
	}


		
}

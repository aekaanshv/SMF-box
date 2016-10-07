#include "krig.h"
double mean_cal(int num_pts, MatD L, VecD y0)
{
//............ CALCULATE MEAN 
        VecD ones(num_pts,1.0),mean_dummy(num_pts);
        chol_sol(num_pts,L,ones,mean_dummy);
        double mean_denom = std::accumulate(mean_dummy.begin(), mean_dummy.end(), 0.0);
        chol_sol(num_pts,L,y0,mean_dummy);
        double mean_num = std::accumulate(mean_dummy.begin(), mean_dummy.end(), 0.0);
        return mean_num/mean_denom;
}

double sigma2_cal(int num_pts, MatD L, VecD y0, double mean)
{

//............ CALCULATE SIGMA^2
	int i;

        VecD y0_minus_mean(num_pts),var_dummy(num_pts);
	for(i=0;i<num_pts;i++)
		y0_minus_mean[i] = y0[i] - mean;

        chol_sol(num_pts,L,y0_minus_mean,var_dummy);
	double total = 0;

	for(i=0;i<num_pts;i++)
		total = total + y0_minus_mean[i]*var_dummy[i];
	return total/num_pts;
}


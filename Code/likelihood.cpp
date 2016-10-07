#include "krig.h"

// DETERMINANT OF A TRIANGULAR MATRIX
// IS SIMPLY THE PRODUCT OF THE DIAGONAL ELEMENTS.


double cal_det_tri(MatD L,int n)
{
	double det=1.0;
	for (int i=0;i<n;i++)
		det *= L[i][i];
	return det;
}

bool bndchk(int dim, VecD points, VecD amin, VecD amax)
{
        for(int i=0; i<dim ; i++)
        {
                if((points[i] <= amin[i]) || ( points[i] >= amax[i]))
                        return false;
        }
        return true;

}

double loglikelihood(model S, double const *theta_d)
{
	
	int num_pts,dim,dim_th;
	VecI th_index; 
	MatD x0,L; 
	VecD y0;
	double MAX = 9E18;
	dim = S.dim;
	dim_th = S.dim_th;
	VecD theta(dim_th);

	for(int i=0; i<dim_th;i++)
		theta[i] = theta_d[i];
	S.theta = theta;
	if(!model2param_wchol(S,num_pts,dim,dim_th,x0,y0,th_index,theta,L))
		return MAX;


//	if(!bndchk(dim,theta,thlow,thupp))
//		return MAX;
	double det = 0;
// DETERMINANT OF FULL MATRIX R IS SIMPLY SQUARE OF DETERMINANT OF L	
	det = cal_det_tri(L,num_pts);
	det = det*det;

// CALCULATE MEAN
  	double mean=0; 
  	mean =  mean_cal(num_pts,L,y0);
	double detpad = 1e-10;
// CALCULATE SIGMA2
   	double sigma2=0;
   	sigma2 = sigma2_cal(num_pts,L,y0,mean);

//	if(sigma2 <= 0 )
//		std::cout<<"Neg sigma hit \n";

// CALCULATE LIKELIHOOD
	double likeli = 0.5*(num_pts*log(sigma2) + num_pts + log(det+detpad));

// ADD PENALTY FOR LARGE THETA :
	double k_penalty = 1.0E1;
	double th_sum = 0.0;
	for( int i=0; i<dim_th ; i++)
		th_sum = th_sum + theta[i];
	likeli = likeli + k_penalty*log(th_sum); 

/*	std::cout<<"MEAN  : "<<mean<<"\n";
        std::cout<<"DET  : "<<det<<"\n";
        std::cout<<"SIGMA2  : "<<sigma2<<"\n";
//        std::cout<<"MEAN  : "<<mean<<"\n";

*/
//	std::cout<<" LIKELIHOOD EVAL : "<<likeli-num_pts/2<<"\n";


	return (likeli-num_pts/2); 
}

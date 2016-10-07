#include "krig.h"
#include <algorithm>


// GENERALIZED PATTERN SEARCH :  N+1 VARIANT
// FOR FINDING A MAXIMUM TO A TOLERANCE

// MADE FOR OPTIMIZING THE MAXIMUM LIKELIHOOD FUNCTION
// ASSUMING IT IS UNIMODAL

/*
bool bndchk(int dim, VecD points, VecD amin, VecD amax)
{
	for(int i=0; i<dim ; i++)
	{
		if((points[i] <= amin[i]) || ( points[i] >= amax[i]))
			return false;
	}
	return true;

}
*/


void vec_copy(VecD x, VecD &y, int n)
{

// 	COPY ENTRIES OF VECTOR X IN TO VECTOR Y
	for (int i=0;i<n;i++)
		y[i]= x[i];
}

void genpattsearch(int dim,VecD amin,VecD amax, std::vector<int> spc, VecD &thmax, double &Jmax, model S, double (*feval)(model))
{

//	std::cout<<" Entered generalized pattern search using MADS \n ";
	int tot_feval = 0;
        int maxpos,max_hist = 4000;
	double MAX_NEG = -9E20;

	VecD delta(dim),th0(dim),J(dim+1);
	double del_sum,tol;
	MatD th_hist,poll_pts;
	VecD J_hist(max_hist),temp(dim);


	int i,thptr=0;

	tol = 0.00000001;


/*	INITIALIZE FIRST POINT - HERE TAKEN AT CENTER OF THE DOMAIN

	for (i=0;i<dim;i++)
		th0[i]= (amax[i]-amin[i])/50 + amin[i];
*/

//	MANUAL INITIALIZATION FROM CALLING FUNCTION
	th0 = thmax;
	S.theta = th0;

//	EVALUATE AT INITIAL POINT
	J_hist[thptr] = (*feval)(S);
	tot_feval++;
	Jmax = J_hist[thptr];
	vec_copy(th0,temp,dim);
	th_hist.push_back(temp);
	thptr++;
//  POLL AROUND th0
	poll_mads(th0,dim,spc,amin,amax,poll_pts);

	del_sum = 0;
	for(i=0;i<dim;i++)
	{
		delta[i] = (amax[i]- amin[i])/(spc[i]-1.0);
		del_sum += delta[i];
	}


//	ENTER THE OPTIMIZATION LOOP

	while(del_sum > dim*tol)
	{
		maxpos = -1;

//	FEVAL AT POLL_POINTS
		for(i=0;i<dim+1;i++)
		{
//	BOUNDARY CHECK
			if(bndchk(dim,poll_pts[i],amin,amax))
			{
				S.theta = poll_pts[i];
				J[i] = (*feval)(S);
				tot_feval++;
			}
			else
				J[i] = MAX_NEG; 
			
//	CHECK IF NEW MAXIMUM IS ACHIEVED
			if(J[i] > Jmax)
			{
				Jmax = J[i];
				maxpos = i;
			}
		}
//	IF NEW MAXIMUM ACHIEVED, UPDATE HISTORY VECTORS
		if(maxpos>=0)
		{
			J_hist[thptr] = Jmax;
			vec_copy(poll_pts[maxpos],temp,dim);
			th_hist.push_back(temp);
			thptr++;
		}
//	OTHERWISE REFINE
		else
		{
			del_sum = 0;
			for(i=0;i<dim;i++)
			{
				spc[i] = (spc[i]-1)*4 + 1;
				delta[i] = delta[i]/4.0;
				del_sum += delta[i];
			}
//			std::cout<<"Del_sum :  "<<del_sum<<std::endl; 
		}
//	POLL ONCE MORE
		poll_pts.clear();
		poll_mads(th_hist[thptr-1],dim,spc,amin,amax,poll_pts);
//		std::cout<<" Total feval till now : "<<tot_feval<<std::endl;
	}


	thptr--;

	Jmax = J_hist[thptr];
	vec_copy(th_hist[thptr],thmax,dim);

//	POST PROCESS
	std::cout<<"Total Likelihood evaluations :	"<<tot_feval<<std::endl;
	std::cout<<"Number of improvements :   "<<thptr<<std::endl;
	std::cout<<" \nPosition of maximum hyperparameter :    "<<std::endl;
	for (i = 0; i < dim; i++)
	{
		std::cout<<thmax[i]<<"   ";
	}
	std::cout<<" \n \n ";
	std::cout<<" \n Final cost maximum :    "<<Jmax<<std::endl;


}

/*DEBUG : TEST MAIN

int main()
{
	int n = 8;
	VecD xmax(n),amin(n),amax(n);
	double Jmax;
	int i,j;
	std::vector<int> spc(n);

	for(i=0;i<n;i++)
	{
    		amin[i] = -35;
		amax[i] = 12;
		spc[i] = 21;
    	}

	std::cout<<"Before opt \n";

	model S; // DEFINE SOMETHING FOR S
	genpattsearch(n,amin,amax,spc,xmax,Jmax,S,loglikelihood);

}*/

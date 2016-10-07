#include "krig.h"

void model2param(model S,int &num_pts, int &dim, int &dim_th, MatD &x0, VecD &y0,VecI &th_index, VecD &theta, MatD &L)
{

        num_pts = S.num_pts;
        dim = S.dim;
	dim_th = S.dim_th;
        x0 = S.x0;
        y0 = S.y0;
	th_index = S.th_index;
        theta = S.theta;
        L = S.L;

        return;
}


bool model2param_wchol(model S,int &num_pts, int &dim, int &dim_th, MatD &x0, VecD &y0,VecI &th_index, VecD &theta, MatD &L)
{

	MatD R;
        num_pts = S.num_pts;
	dim_th = S.dim_th;
        dim = S.dim;
        x0 = S.x0;
        y0 = S.y0;
	th_index = S.th_index;
        theta = S.theta;
        VecD zeros(num_pts,0);
        L.clear();
        int i;
// GENERATE COVARIANCE MATRIX

	for(i=0;i<num_pts;i++)
        {
               	VecD V;
                expcov(num_pts,dim,dim_th,x0,x0[i],th_index,theta,V);
       	        R.push_back(V);
                L.push_back(zeros);
       	}
// DECOMPOSE COVARIANCE MATRIX
       	 if(!Cholesky(num_pts,R,L))
        {
//                std::cout<<"Cholesky failed - negative factors in square root"<<std::endl;
        	return false;
	}
	S.L = L;
	return true;
}

bool param2model(model &S, int num_pts,int dim, int dim_th, MatD x0, VecD y0,VecI th_index, VecD theta)
{
        VecD zeros(num_pts,0);
        S.num_pts = num_pts;
        S.dim = dim;
	S.dim_th = dim_th;
        S.x0 = x0;
        S.y0 = y0;
	S.th_index = th_index;
        S.theta = theta;
	S.L.clear();
        int i;
	MatD R,L;
// GENERATE COVARIANCE MATRIX
	R.clear();
	L.clear();
        for(i=0;i<num_pts;i++)
        {
                VecD V;
                expcov(num_pts,dim,dim_th,x0,x0[i],th_index,theta,V);
                R.push_back(V);
                L.push_back(zeros);
        }
// DECOMPOSE COVARIANCE MATRIX
         if(!Cholesky(num_pts,R,L))
        {
//                std::cout<<"Cholesky failed - negative factors in square root"<<std::endl;
                return false;
        }
        S.L = L;
	return true;

}

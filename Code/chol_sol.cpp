#include "krig.h"

// SOLVE AX = B
// WHERE A HAS BEEN PRE-FACTORIZED INTO L*L'

void fwd_sub(int num_pts,MatD L, VecD B, VecD &X)
{

// SOLVE Lx = B

	double sum;
	for(int i=0;i<num_pts;i++)
	{
		sum = 0.0;
		for(int j=0;j<i;j++)
			sum = sum + L[i][j]*X[j];
		X[i] = ((B[i] - sum)/L[i][i]);
	}
}


void back_sub(int num_pts,MatD L, VecD B, VecD &X)
{

// SOLVE L'x = B

        double sum;
        for(int i=num_pts-1;i>=0;i--)
        {
                sum = 0.0;
                for(int j=i+1;j<num_pts;j++)
                        sum = sum + L[j][i]*X[j];
                X[i] = ((B[i] - sum)/L[i][i]);
        }
}

void chol_sol(int num_pts,MatD L, VecD B, VecD &X)
{
	VecD B1(num_pts);
	fwd_sub(num_pts,L,B,B1);
	back_sub(num_pts,L,B1,X);
}




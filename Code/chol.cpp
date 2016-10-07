#include "krig.h" 


// Cholesky requires the matrix to be symmetric positive-definite
// Returns D = Lower triangular dense matrix. Can be optimized to conserve memory by converting D into a single vector.

bool Cholesky(int d,MatD R,MatD &L)
{
   double sum,fact;
   int i,j,k,p; // Dummy indices
   for(k=0;k<d;k++)
	for(p =0;p<d;p++)
		L[k][p] = 0.0;
  
// CHECK FOR SIMILAR ROWS IN R. IF FOUND, NUGGET THEM
// Update this part to the numeric accumulator 

   double norm; // 2-Norm of difference between rows
   double tol = 1e-10;
   double nugget = 1e-8;
   std::vector<int> nugget_these_rows(d,0);
   for(i=0;i<d;i++)   
   {
	VecD V;		
	for(k=i+1;k<d;k++) 
	{
		norm=0;
		for(j=0;j<d;j++) 
			V.push_back(R[k][j] - R[i][j]);
		for(j=0;j<d;j++)
			norm += V[j]*V[j];
		if(norm < tol*tol)
			nugget_these_rows[k] = 1;
	}
   }
   for(i=0;i<d;i++)
   {
	if(nugget_these_rows[i]==1)
	{
		R[i][i] = R[i][i] + nugget;
	}
   }
	
// CALCULATE CHOLESKY DECOMPOSITION
   for(k=0;k<d;++k)
   {
      sum=0.0;
      for(p=0;p<k;++p)
         sum+=L[k][p]*L[k][p];
      fact = R[k][k]-sum;
      if(fact > tol*tol )
         L[k][k]=sqrt(fact);
      else
         return false;      
      for(i=k+1;i<d;++i)
      {
         sum=0.;
         for(p=0;p<k;++p)
            sum+=L[i][p]*L[k][p];
         L[i][k]=(R[i][k]-sum)/L[k][k];
      }
   }
   return true;

}

void expcov(int num_pts,int dim,int dim_th,MatD x,VecD x0,VecI th_index,VecD theta,VecD &covx)
{
	const double e = 2.7182818284;
	double sum;
	double factor;
//
	if(dim < 8 )
		factor = 1;
	else
		factor = 1.1/dim + 10.1/(dim*dim);
//
	for(int i=0;i<num_pts;i++)
	{
		sum = 0.0;
		for(int j=0;j<dim;j++)
//			sum = sum + (-theta[j])*(x[i][j]-x0[j])*(x[i][j]-x0[j]);
			sum = sum + factor*(-theta[th_index[j]])*(x[i][j]-x0[j])*(x[i][j]-x0[j]);
		covx.push_back(pow(e,sum));
	}

//	std::cout<<covx.size()<<std::endl;

}

//	DEBUG : MAIN TO TEST CHOLESKY + EXPONENTIAL COVARIANCE

/*

int main()
{

	int i,j,num_pts=4,dim=3;
	MatD x,R,L,inv;

//	INITIALIZATION TO ZERO
	VecD zeros(num_pts,0),theta(dim,2);
	double t;
	for(i=0;i<num_pts;i++)
	{
		VecD temp;
		for(j=0;j<dim;j++)
			temp.push_back(0.001*i+0.0001*j);
		x.push_back(temp);
	}

	std::cout<<"Before exp"<<std::endl;
	for(i=0;i<num_pts;i++)
	{
		VecD V;
		expcov(num_pts,dim,x,x[i],theta,V);
		R.push_back(V);
		L.push_back(zeros);
		inv.push_back(zeros);
	}


	if (Cholesky(num_pts,R,L))
	{
		std::cout<<"L matrix : "<<std::endl;
		for(i=0;i<num_pts;i++)
		{
			for(j=0;j<num_pts;j++)
				std::cout<<L[i][j]<<" \t";
			std::cout<<std::endl;
		}	

//	CALC L'LR. THIS SHOULD BE IDENTITY

		for(i=0;i<num_pts;i++)
		{
			for(j=0;j<num_pts;j++)
			{
				t=0;
				for(int k=0;k<num_pts;k++)
				{
//					for(int p=0;p<n;p++)
						t+=L[i][k]*L[j][k];
				}
				inv[i][j] = R[i][j] - t;
			}
		}

                std::cout<<"Difference between LL' and R : "<<std::endl;
                for(i=0;i<num_pts;i++)
                {
                        for(j=0;j<num_pts;j++)
                                std::cout<<inv[i][j]<<" \t";
                        std::cout<<std::endl;
                }

	}
	else
		std::cout<<"Cholesky failed - negative factors in square root"<<std::endl;
	

}
*/ 


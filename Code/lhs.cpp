#include "krig.h"

double lhs_merit(int dim, int n_div, MatI perm_mat)  // dim x ndiv points
{
	int p = 50;  // phi_p parameters 
	int t = 1;

	double phi_p = 0.0;
	double temp;
	for(int i = 0; i < n_div; i++)
	{
		for(int j = i+1 ; j < n_div ; j++)
		{
			double sumk = 0.0;
			for(int k = 0; k < dim ; k++)
			{	
				temp = pow(fabs((double)(perm_mat[k][i]-perm_mat[k][j])),t);
				sumk += temp;
			}
			sumk = pow(sumk,1.0/t);

			phi_p += pow(sumk,-p);
		}
	}

	phi_p = pow(phi_p,1.0/p);

	return phi_p;

}


void lhs(int dim, int n_div, settings &set0 , bool imp_flag, MatD &x,VecI xlog_index)
{
	int i,j;
//	int n_div = N_DIV1;
	VecI perm(n_div);
	MatI perm_mat;

	for(i=0;i<dim;i++)
	{
		for(j=0;j<n_div;j++)
			perm[j]=j;
		std::random_shuffle(perm.begin(),perm.end(),randn);
		perm_mat.push_back(perm);
	}

// COLUMNS OF PERM_MAT HAVE OUR REQUIRED LHS INDICES

// ---- CODE FOR TESTING INTIAL LHS DESIGN : 

	MatD xtemp;
        for (i = 0; i < n_div; i++)
        {
		
                VecD temp(dim);

                for (j = 0; j < dim; j++)
                {
                        temp[j] = set0.xlow[j] + ((double)(perm_mat[j][i])*(set0.xupp[j]-set0.xlow[j]))/((double)n_div);

                        temp[j] += ((double)(1+randn(set0.spc[j]-2))*(set0.xupp[j]-set0.xlow[j]))/((double)n_div*set0.spc[j]);
                }
                push_to_grid(temp,set0.xlow,set0.xupp,set0.spc);
                xtemp.push_back(temp);
        }
	print_matrix_to_txt_file("x_init_lhs.dat",xtemp,n_div,dim,xlog_index);

// ---- END CODE FOR TSTING INITIAL LHS DESIGN

	if(imp_flag)
	{
		int niter = 10000 + 1000*dim;

		// --- IMPROVE LATIN HYPERCUBIC SAMPLING - SIMULATED ANNEALING
		// --- RANDOMIZE ALONG A PARTICULAR DIRECTION 
		// --- CHECK FOR IMPROVEMENT USING PHI_P MINIMIZATION

		double temperature = 1.0;

		double merit_original = lhs_merit(dim, n_div, perm_mat);
		MatI perm_mat_new = perm_mat;
		MatI perm_mat_best = perm_mat; 
		double merit_best = merit_original; 


		VecD dim_to_improve(dim,1.0);  // -- VECTOR CONTAINING WEIGHTS FOR DIMENSIONS WHERE WE SHOULD FOCUS IMPROVEMENT

		double improve_factor = 1.0 - (1.0/pow(2,dim));


		std::cout<<"Value of original phi_p = "<<merit_original<<std::endl;
	
		for(int iter = 0; iter< niter ; iter++)
		{
			temperature = 1.0 -pow((double)iter/niter,2);
			double sum = 0;

			for(j = 0; j < dim; j++)
				sum = dim_to_improve[j];
			double rand_var = randd(0.0,sum);

			int index = 0;
			while (rand_var >= 0.0)
			{
				rand_var -= dim_to_improve[index]; 
				index++ ;
			}
			index-- ;  // INDEX FOR SHUFFLING 

			// -- SHUFFLE A ROW OF PERM_MAT AND CHECK IF THERE IS IMPROVEMENT

			VecI perm(n_div);
			for(j = 0; j<n_div;j++)
				perm[j] = perm_mat_new[index][j];

			std::random_shuffle(perm.begin(),perm.end(),randn);

			for(j = 0; j<n_div;j++)
				perm_mat_new[index][j] = perm[j] ;

			double merit_new = lhs_merit(dim, n_div, perm_mat_new);

			double ap = exp((merit_new - merit_original)/temperature);

			if(merit_new < merit_original) // BETTER DESIGN FOUND
			{
				perm_mat = perm_mat_new ; 
				merit_original = merit_new ;
				if(merit_new < merit_best) 
				{
                                	perm_mat_best = perm_mat_new ;
                                	merit_best = merit_new ;
				}
				for(j = 0 ; j < dim ; j++)  // RESET IMPROVEMENT VECTOR
					dim_to_improve[j] = 1.0;
			}
			else
			{
				if(ap > randd(0.0,1.0))
				{
					perm_mat = perm_mat_new ;
	                                merit_original = merit_new ;
				}
				dim_to_improve[index] = dim_to_improve[index]*improve_factor;				
			}

			if(iter%1000 == 0)
				std::cout<<"Value of phi_p = "<<merit_best<<" at iter = "<<iter<<std::endl;
		}

		perm_mat = perm_mat_best;

	}

	for (i = 0; i < n_div; i++)
	{
		VecD temp(dim);

                for (j = 0; j < dim; j++)
                {
			temp[j] = set0.xlow[j] + ((double)(perm_mat[j][i])*(set0.xupp[j]-set0.xlow[j]))/((double)n_div);

			temp[j] += ((double)(1+randn(set0.spc[j]-2))*(set0.xupp[j]-set0.xlow[j]))/((double)n_div*set0.spc[j]);
		}
		push_to_grid(temp,set0.xlow,set0.xupp,set0.spc);
		x.push_back(temp);

		// Increment number of lhs points to compute

		set0.npts_lhs += 1;

	}


}

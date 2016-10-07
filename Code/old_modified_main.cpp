#include "krig.h"
#define DIM_1 5   // Dimension of problem 
#define DIM_2 2   // Dimension of hyperparamter space
#define N_DIV1 21 


void set_th_index(VecI th_index)
{

//	MODIFY IN ACCORDANCE TO TH_INDEX THAT YOU WANT	
	for(int i=0;i<DIM_1;i++)
		th_index[i] = 0;

	th_index[2] = 1;
	th_index[3] = 1;
	th_index[4] = 1;

}

void set_th_bounds(VecD thlow,VecD thupp)
{
//	MODIFY IN ACCORDANCE TO TH_BOUNDS THAT YOU WANT
//  ALSO KEEP TRACK OF MAPPING BETWEEN TH_INDEX TO ACTUAL THETA	
	for(int i=0;i<DIM_2;i++)
	{
		thlow[i] = 1e-2;
		thupp[i] = 2e1;		
	}
}

void set_grid(VecD xlow,VecD xupp, VecI spc)
{
//	MODIFY IN ACCORDANCE TO THE CONSTRAINTS AND GRID SPACING THAT YOU WANT

	for(int i=0;i<DIM_1;i++)
	{
		xlow[i] = -5;
		xupp[i] = 10;		
		spc[i] = 117;
	}
}

void lhs(int dim, VecD xlow,VecD xupp, VecI spc, VecD x0, VecD x)
{
	int i,j;
	int n_div = N_DIV1;
	VecI perm(n_div);
	std::vector< std::vector<int> > perm_mat;

	for(i=0;i<dim;i++)
	{
		for(j=0;j<n_div;j++)
			perm[j]=j;
		std::random_shuffle(perm.begin(),perm.end(),randn);
		perm_mat.push_back(perm);
	}

// COLUMNS OF PERM_MAT HAVE OUR REQUIRED LHS INDICES


	for (i = 0; i < n_div; i++)
	{
		VecD temp(dim);

                for (j = 0; j < dim; j++)
                {
			temp[j] = xlow[j] + ((double)(perm_mat[j][i])*(xupp[j]-xlow[j]))/((double)n_div);

			temp[j] += ((double)(1+randn(spc[j]-2))*(xupp[j]-xlow[j]))/((double)n_div*spc[j]);
		}
		push_to_grid(temp,xlow,xupp,spc);
		x0.push_back(temp);
		x.push_back(temp);
	}

}


double blackbox(VecD x)
{
// ENTER BLACK BOX FUNCTION HERE
	int n = x.size();
	double tot = 0;

 /*  QUADRATIC + BUMP
	for(int i=0;i<n;i++)
		tot += -1.0*(x[i] - x[i]*x[i]/2) - exp(-20*(x[i]-1)*(x[i]-1)); ;  
//	tot = tot*n;

*/
// ROSENBROCK FUNCTION 

	for(int i=0;i<n-1;i++)
	{
		tot += pow((1-x[i]),2) + 100*pow((x[i+1] -x[i]*x[i]),2);	
	}
//
	return tot;
}


bool newpoint(int num_pts, int dim, MatD x0, VecD x_cand)
{
       	double x_dist = 0.0;
        double x_dist_tol = 1e-6;
	for(int i =0;i<num_pts;i++)
        {
		x_dist = 0;
		for(int j=0;j<dim;j++)
               		x_dist += pow(x0[i][j] - x_cand[j],2);
		x_dist = pow(x_dist,0.5);
        	if(x_dist < x_dist_tol)
			return false;
	}
	return true;
		
}



int main()
{

//	SEED RANDOM GENERATOR FIRST
	std::srand(static_cast <unsigned> (time(0)));

    double Jmin=1e10,Jnew,Jkrig;
	int i,j,num_pts=0;
	MatD x;
	VecD J;

	int dim = DIM_1;
	int dim_th = DIM_2;

	VecD xmin(dim),xkrig(dim);

//	THETA INDEX VECTOR
	VecI th_index(dim);
//  SET BOUNDS FOR THETA
	VecD thlow(dim_th);
	VecD thupp(dim_th);
	VecD theta(dim_th);

	set_th_index(th_index);
	set_th_bounds(thlow,thupp);

//  INITIALIZE THETA

	for(i=0;i<dim_th;i++)
		theta[i]= sqrt(thupp[i]*thlow[i]);


	FILE *fp;
	fp = fopen("BestJHist.dat","w");


//=============== DEFINE COUNTERS FOR SMF =====================

	int  n_poll=0, n_search=0, n_refine=0;

//=============== DEFINE BOX CONSTRAINTS ======================

	VecD xupp(dim);
	VecD xlow(dim);
	VecI spc(dim); 

	set_grid(xlow,xupp,spc);

//=============== SMF PARAMETERS ==============================

	double conv_tol = 0.000005;
	double x_del_sum = 0;
	VecD x_delta(dim);

	for(i=0;i<dim;i++)
	{
		x_delta[i] = (xupp[i] - xlow[i])/(spc[i]-1);
		x_del_sum += x_delta[i];
	}	

//=============== READ CURRENT STATE IN SMF LOOP =============

	int smf_phase;
	read_phase("smfphase.dat",smf_phase);

	if(phase == 1)		// LHS
	{
		// Initialize some points where function will be evaluated
// Use LHS sampling 

		lhs(dim,xlow,xupp,spc,x);

/*		PRINT LHS POINTS 

		std::cout<<" \n \n Finally after lhs : \n \n"; 
	
		for(i=0;i<n_div;i++)
		{
			for(j=0;j<dim;j++)
				std::cout<<x[i][j]<<"  ";
			std::cout<<std::endl;
		}
*/


/*		FINDING XMIN

	for(i=0;i<n_div;i++)
	{
		J.push_back(blackbox(x[i]));
		n_feval++;
		if(J[i] < Jmin)
		{
			Jmin = J[i];
			xmin = x[i];
		}
		write_hist(fp,xmin,Jmin);	
	}
	num_pts = J.size();
	std::cout<<" Number of LHS feval : "<<n_feval<<std::endl;
	std::cout<<" Jmin from LHS : 	   "<<Jmin<<std::endl;

*/
// ============== EVALUATE COST FUNCTION AT THOSE POINTS ======

		num_pts = N_DIV1;
		print_matrix_to_txt_file("pts_to_eval.dat",x,num_pts,dim);
		write_phase("smfphase.dat","SEARCH");
		return;

	}

	else if(phase == 2)		// SEARCH PHASE
	{
		int  num_pts_y;
		import_matrix_from_txt_file("xhist.dat",x,num_pts,dim);
		import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
		model S;
		double likeli_max;
		if(num_pts == num_pts_y)
			std::cout<<" \nNumber of function evaluations till now :   "<<num_pts<<"\n \n";
		else
		{
			std::cout<<" \nInput error : number of points in xhist dont match Jhist \n";
			return;
		}

// ============== GENERATE SURROGATE FOR THE POINTS ===========


		param2model(S,num_pts,dim,dim_th,x,J,th_index,theta);
		int nbrestarts = 0;
		double incpopsize = 2.0;
		double *theta_d;
		double tol_th = 1e-3;
		VecD theta_stddev(dim_th,0.1);
		std::cout<<" HYPER PARAMETER OPTIMIZATION : \n ";
		theta_d =cm_optimize(dim_th, S, loglikelihood, thlow, thupp, nbrestarts, incpopsize,theta,theta_stddev,tol_th);  
		std::cout<<" THETA OBTAINED : \n";

		for(i=0;i<dim_th;i++)
		{
			theta[i] = theta_d[i];
			std::cout<<theta[i]<<"  ";
		}
		std::cout<<"\n \n";

/*		DEBUG : CHECK CHOLESKY OF OPTIMAL

		if(!param2model(S,num_pts,dim,dim_th,x,J,th_index,theta))
        {
                std::cout<<" PARAMETER TO MODEL - FAILED: CHOLESKY. \n \n";
                return 0;
        }

*/
		double mse = msecheck(S);
		std::cout<<"\nMSE of new model=  "<<mse<<std::endl;

/*		DEBUG : CHECK IF MSE IS NORMAL OR NOT

		if(!std::isnormal(mse))
		{
			std::cout<<" DEBUG MODE : MSE TURNED ABNORMAL -.- \n ";
			std::cout<<" OUTPUT MODEL : (X0)  feval : Y0 \n ";

			std::cout<<" ----------------------- \n ";
			for(i=0;i<num_pts;i++)
			{
				for(j=0;j<dim;j++)
					std::cout<<S.x0[i][j]<<"   ";
				std::cout<<"    feval : "<<S.y0[i]<<" \n";
			}
			std::cout<<" \n THETA :\n";
			for(i=0;i<dim;i++)
				std::cout<<S.theta[th_index[i]]<<"   ";
			std::cout<<" \n";
			std::cout<<"Numpts  :  "<<S.num_pts<<"\n";
			std::cout<<"Dim  :  "<<S.dim<<"\n";
			std::cout<<"Loglikelihood  : "<<loglikelihood(S,theta_d)<<"\n";

			return 0;
 
		}

*/

//============== CALCULATE NEW CANDIDATE MINIMUM ==============

        nbrestarts = 0;
        incpopsize = 2.0;
        double *xkrig_d;
		n_search++;
		double tolx = 1e-7;
		std::cout<<"SEARCH STEP \n";
		VecD x_stddev(dim,0.7);
        xkrig_d =cm_optimize(dim, S, pred_fn, xlow, xupp, nbrestarts, incpopsize, xmin, x_stddev, tolx);


        for(i=0;i<dim;i++)
               xkrig[i] = xkrig_d[i];


		push_to_grid(xkrig,xlow,xupp,spc);

//=============== CHECK IF NEW POINT IS ACHIEVED ==============

        bool search,poll;
        bool newmin = true;
        num_pts = J.size();
        newmin = newpoint(num_pts,dim,x,xkrig);

//=============== CHECK IF NEW MINIMUM IS ACHIEVED ============
		if(newmin == true)
		{
			x.push_back(xkrig);
			print_matrix_to_txt_file("pts_to_eval.dat",xkrig,1,dim);
			num_pts++;			
			print_matrix_to_txt_file("xhist.dat",x,num_pts,dim);
			return;
		}
		else
			write_phase("smfphase.dat","POLL");

	}

	else if(phase ==3)
	{
			int  num_pts_y;
			import_matrix_from_txt_file("xhist.dat",x,num_pts,dim);
			import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
			if(num_pts == num_pts_y)
				std::cout<<" \nNumber of function evaluations till now :   "<<num_pts<<"\n \n";
			else
			{
				std::cout<<" \nInput error : number of points in xhist dont match Jhist \n";
				return;
			}

			MatD poll_points,poll_to_eval;
			//-----------insert code for xmin-------------
        	poll_mads(xmin,dim,spc,xlow,xupp,poll_points);
        	for(i=0;i<dim+1;i++)
			{
				if(newpoint(num_pts,dim,x,poll_points[i]))
				{
					poll_to_eval.push_back(poll_points[i]);
				}
			}


	}

// ============== BEGIN SMF LOOP ==============================



	while(x_del_sum > dim*conv_tol)
	{






		if((Jkrig < Jmin)&&(newmin == true))
		{
//=============== SUCCESSFUL SEARCH : UPDATE POINTS ===========
			std::cout<<"Search successful \n";
			Jmin = Jkrig;
			xmin = xkrig;
	                search = true;
		}
		else 
			search = false;
	
		write_hist(fp,xmin,Jmin);





    if(search == false)
	{
//=============== FAILED SEARCH : POLL AROUND MINIMUM =========
        	std::cout<<"Search failed \n";
        	search = false;
        	poll_mads(xmin,dim,spc,xlow,xupp,poll_points);

//=============== CHECK IF POLL GIVES NEW MINIMUM =============
        	n_poll++;
        	poll = false;
        	double J_temp;

 		for(i=0;i<dim+1;i++)
		{
			if(newpoint(num_pts,dim,x,poll_points[i]))
			{
				J_temp = blackbox(poll_points[i]);
				n_feval++;
			
//=============== SUCCESSFUL POLL : UPDATE POINTS =============

				if(J_temp < Jmin)
				{
					std::cout<<" Poll step improvement : Jmin old = "<<Jmin<<"   and new  : "<<J_temp<<std::endl;
					Jmin = J_temp;
					poll = true;
					xmin = poll_points[i];
				}
				J.push_back(J_temp);
				x.push_back(poll_points[i]);
				write_hist(fp,xmin,Jmin);
			}
		}       	

//=============== FAILED POLL : REFINE ========================
		if(poll == false)
		{
			std::cout<<"Poll failed \n";
			n_refine++;
			x_del_sum = 0;
			for(i=0;i<dim;i++)
			{
                		spc[i] = (spc[i]-1)*4 + 1;
                		x_delta[i] = x_delta[i]/4.0;
                		x_del_sum += x_delta[i];
            }
		}
		else
			std::cout<<"Poll worked. \n \n";

        }
	std::cout<<" X_DEL_SUM =  "<<x_del_sum<<std::endl;
//=============== END SMF LOOP ================================

	} // END OF SMF

	fclose(fp);
//=============== POST-PROCESS ================================
	std::cout<<"Size of J 	     :			"<<J.size()<<std::endl;
	std::cout<<"Minimum achieved :			"<<Jmin<<std::endl;
	std::cout<<"Position of minima :		"<<std::endl;
	for(i=0;i<dim;i++)
		std::cout<<xmin[i]<<"  ";
	std::cout<<"\n \n \n";
	std::cout<<"Number of function evaluations :  "<<n_feval<<std::endl;
	std::cout<<"Number of search steps :  	"<<n_search<<std::endl;
	std::cout<<"Number of poll steps :  	"<<n_poll<<std::endl;
	std::cout<<"Number of refine steps : 	"<<n_refine<<std::endl;
}

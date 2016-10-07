#include "krig.h"
#define DIM_1 6   // Dimension of problem 
#define DIM_2 DIM_1 + 0   // Dimension of hyperparamter space
#define N_DIV1 3*DIM_1 + 1 

void set_th_index(VecI &th_index)
{

//	MODIFY IN ACCORDANCE TO TH_INDEX THAT YOU WANT	
	for(int i=0;i<DIM_1;i++)
		th_index[i] = i;
//	th_index[1] = 1;
/*  abhay
	th_index[1] = 1;
	th_index[3] = 1;
	th_index[5] = 1;
	th_index[7] = 2;
*/
}

void set_xlog(VecI &xlog_index)
{
	int dim = xlog_index.size();
//	xlog_index[1] = 1;

	for(int i =0;i<dim;i++)
		xlog_index[i]=0;	
/*	xlog_index[0] =1;
	xlog_index[2] =1;
	xlog_index[4] =1;
	xlog_index[6] =1;
*/
}

void set_th_bounds(VecD &thlow,VecD &thupp)
{
//	MODIFY IN ACCORDANCE TO TH_BOUNDS THAT YOU WANT
//  ALSO KEEP TRACK OF MAPPING BETWEEN TH_INDEX TO ACTUAL THETA	
	for(int i=0;i<DIM_2;i++)
	{
		thlow[i] = 1e-4;
		thupp[i] = 200;		
	}

/*	thlow[0] = 1E-9;
	thlow[2] = 1E-9;
	thlow[4] = 1E-9;
	thlow[6] = 1E-9;

	thupp[0] = 1E-9;
        thupp[2] = 1E-9;
        thupp[4] = 1E-9;
        thupp[6] = 1E-9;
*/	 


//	std::cout<<"Theta bounds set \n";
}

// Old bounds read routine
// void set_grid(VecD &xlow,VecD &xupp, VecI &spc, VecI xlog_index)
// {
// //	MODIFY IN ACCORDANCE TO THE CONSTRAINTS AND GRID SPACING THAT YOU WANT
// 	int dim = DIM_1;

	
// 	VecD spc_d,xlow_d,xupp_d;
// 	int dim_chk;
// 	import_vector_from_txt_file("xlow.dat",xlow_d,dim_chk);

//         if(dim_chk != dim)
//                 std::cout<<"INCORRECT XLOW FILE : WRONG DIMENSION!! \n";
// 	xlow = xlow_d;

// 	import_vector_from_txt_file("xupp.dat",xupp_d,dim_chk);
//         if(dim_chk != dim)
//                 std::cout<<"INCORRECT XUPP FILE : WRONG DIMENSION!! \n";
// 	xupp = xupp_d;

// 	import_vector_from_txt_file("spc.dat",spc_d,dim_chk);


// 	if(dim_chk != dim)
// 		std::cout<<"INCORRECT SPC FILE : WRONG DIMENSION!! \n";
// 	for(int i=0;i<dim;i++)
// 	{
// 		spc[i] = (int)spc_d[i];
// 	}

// 	for(int i=0;i<dim;i++)
// 	{
// 		if(xlog_index[i]==1)
// 		{
// 			xlow[i] = log10(xlow[i]);
// 			xupp[i] = log10(xupp[i]);
// 		}
// 	}

// }

void find_xmin(MatD x,VecD J, VecD &xmin, double &Jmin)	// Find set of parameters for current incumbent solution
{
	int num_pts = J.size();
	int min_pos = 0;
	Jmin = 1E10;
	for(int i=0;i<num_pts;i++)
	{
		if(Jmin > J[i])
		{
			Jmin = J[i];
			min_pos = i;
		}
	}

	xmin = x[min_pos];
}

bool newpoint(int num_pts, int dim, MatD x0, VecD x_cand)
{
        double x_dist = 0.0;
        double x_dist_tol = 1e-7;
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

void terminate_opt(int terminate)
{
        VecD term(1);
        term[0] = (double)terminate;

        print_vector_to_txt_file("stopfile",term,1);

}

bool succ_search(int dim,settings &set0,VecI xlog_index)
{
	MatD x,xnew;
	VecD J,Jnew;
	int num_pts,num_pts_y,num_pts_new;
	bool is_success;
	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
	import_matrix_from_txt_file("xnew.dat",xnew,num_pts_new,dim,xlog_index);
	import_vector_from_txt_file("Jnew.dat",Jnew,num_pts_new);
	VecD xmin,xmin_new;
	double Jmin,Jmin_new;
	find_xmin(x,J,xmin,Jmin);
	find_xmin(xnew,Jnew,xmin_new,Jmin_new);


	if(Jmin_new < Jmin)
	{
		is_success = true;
		FILE *fp;
		fp = fopen("search_stats","a");
		fprintf(fp,"%.5f \n",Jmin-Jmin_new);
		fclose(fp);

		set0.succ_search += 1 ;

	}
	else
		is_success = false;

	for (int i = 0; i < num_pts_new; i++)
	{
		x.push_back(xnew[i]);
		J.push_back(Jnew[i]);
	}


	num_pts += num_pts_new;
	num_pts_y += num_pts_new;



	print_matrix_to_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	print_vector_to_txt_file("Jhist.dat",J,num_pts_y);
	return is_success;


}

bool succ_poll(int dim, settings &set0, VecI xlog_index)
{
    MatD x,xnew;
    VecD J,Jnew;
    int num_pts,num_pts_y,num_pts_new;
	bool is_success;
	
	// Read history and candidate data

	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
	import_matrix_from_txt_file("xnew.dat",xnew,num_pts_new,dim,xlog_index);
	import_vector_from_txt_file("Jnew.dat",Jnew,num_pts_new);
	
	// Check if new minimum is achieved

	VecD xmin,xmin_new;
	double Jmin,Jmin_new;
	find_xmin(x,J,xmin,Jmin);
	find_xmin(xnew,Jnew,xmin_new,Jmin_new);


	if(Jmin_new < Jmin)
	{
		is_success = true;
		FILE *fp;
        fp = fopen("poll_stats","a");
        fprintf(fp,"%.5f \n",Jmin-Jmin_new);
        fclose(fp);

		set0.succ_poll += 1 ;


	}
	else
		is_success = false;

	for (int i = 0; i < num_pts_new; i++)
	{
		x.push_back(xnew[i]);
		J.push_back(Jnew[i]);
	}

	// Update counters

	num_pts += num_pts_new;
	num_pts_y += num_pts_new;

	print_matrix_to_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	print_vector_to_txt_file("Jhist.dat",J,num_pts_y);
	return is_success;
}

bool search(int dim, int dim_th, VecD thlow, VecD thupp,VecI th_index, settings &set0, VecI xlog_index)
{

	MatD x,xnew;
	VecD J,Jnew,theta;

	VecD xkrig(dim);

	int  num_pts,num_pts_new,num_pts_y,dim_th_chk;
	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);

//	import_vector_from_txt_file("thopt.dat",theta,dim_th_chk);	

	theta = set0.theta;


// 	if(dim_th_chk != dim_th)
// 	 	std::cout<<"ERROR IN SETTINGS - THETA HAS WRONG DIMENSIONS \n"<<dim_th_chk<<std::endl; 

	model S;
	double likeli_max;
	if(num_pts == num_pts_y)
		std::cout<<" \nNumber of function evaluations till now :   "<<num_pts<<"\n \n";
	else
	{
		std::cout<<" \nInput error : number of points in xhist dont match Jhist \n";
		return false;
	}

	VecD xmin(dim);

	double Jmin;
	find_xmin(x,J,xmin,Jmin);

	
// ============== GENERATE SURROGATE FOR THE POINTS ===========


	while(!param2model(S,num_pts,dim,dim_th,x,J,th_index,theta))
	{
		std::cout<<" PARAMETER TO MODEL - FAILED: CHOLESKY. \n \n";
		for(int p = 0 ; p < dim_th ; p++)
	    {
		
			if(theta[p] > thupp[p])
			{
				terminate_opt(-1);	
	       		return false;
			}
			else
				theta[p] = theta[p]*pow(2,(1.0/(dim-1)));
	    }
    }


	int nbrestarts = 0;
	double incpopsize = 4.0;
	double *theta_d;
	double tol_th = 1e-4 + 0.00001*num_pts;
	VecD theta_stddev(dim_th,0.1);
	std::cout<<" HYPER PARAMETER OPTIMIZATION : \n ";
	theta_d =cm_optimize(dim_th, S, loglikelihood, thlow, thupp, nbrestarts, incpopsize,theta,theta_stddev,tol_th);  
	std::cout<<" THETA OBTAINED : \n";

	for(int i=0;i<dim_th;i++)
	{
		theta[i] = theta_d[i];
		std::cout<<theta[i]<<"  ";
	}
	std::cout<<"\n \n";

        if(!param2model(S,num_pts,dim,dim_th,x,J,th_index,theta))
        {
                std::cout<<" PARAMETER TO MODEL - FAILED: CHOLESKY. \n \n";
                return false;
        }


	double mse = msecheck(S);
	std::cout<<"\nMSE of new model=  "<<mse<<std::endl;

	set0.theta = theta ; 

//============== CALCULATE NEW CANDIDATE MINIMUM ==============

	// Multi-start search -- Do multiple searches by large scale perturbations to the current best
	
    nbrestarts = 0;
    incpopsize = 4.0;
    double *xkrig_d;
	double tolx = 1e-7;
	std::cout<<"SEARCH STEP \n";

	int n_test_pts = 8;
	double currbest = 1E10 ;

	VecD xmin_pert = xmin ; 

	double stdfact = 0.25 ;

	VecD stdvec(dim);

	for(int k=0; k< dim ;k++)
	{
		stdvec[k] = stdfact*(set0.xupp[k] - set0.xlow[k]);	
	}

	for(int k =0 ; k< n_test_pts ; k++)
	{

		VecD x_stddev = stdvec ;

		if(k>0)
		{
			
			for(int i=0;i<dim;i++)
			{
				xmin_pert[i] = xmin[i] +randd(-2*stdvec[k],2*stdvec[k]) ;  
			}
		}

		VecD tmpx(dim);
		VecD tmpout(1); 
   		xkrig_d =cm_optimize(dim, S, pred_fn, set0.xlow, set0.xupp, nbrestarts, incpopsize, xmin, x_stddev, tolx);

 	  	for(int i=0;i<dim;i++)
 	   		tmpx[i] = xkrig_d[i];

 	   	MatD xmattmp ; 
 	   	xmattmp.push_back(tmpx) ;

 	   	predictor(S,1, xmattmp,tmpout);

 	   	if(tmpout[0]- currbest< -0.001*abs(currbest))
 	   	{
 	   		xkrig = tmpx;
 	   		currbest = tmpout[0]; 
 	   		std::cout<<"Multipoint CMAES : Improved candidate point found \n";
 	   	}
   	}

    push_to_grid(xkrig,set0.xlow,set0.xupp,set0.spc);


//=============== CHECK IF NEW MINIMUM IS ACHIEVED ============
	set0.nstage_search += 1 ; 

	if(newpoint(num_pts,dim,x,xkrig) == true)
	{
		MatD temp;
		x.push_back(xkrig);
		temp.push_back(xkrig);
		print_matrix_to_txt_file("xnew.dat",temp,1,dim,xlog_index);
		std::cout<<"New search point: \n"; 
		for(int i=0; i<dim;i++)
		{
			std::cout<<xkrig[i]<<"  "; 
		}
		std::cout<<"\n \n";

		set0.npts_search +=1 ; 

		return true;
	}
	else
	{	return false;}

}

void poll(int dim, settings &set0, VecI xlog_index)
{
	MatD x;
	VecD J;
	int  num_pts_y,num_pts;
	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
	if(num_pts == num_pts_y)
		std::cout<<" \nNumber of function evaluations till now :   "<<num_pts<<"\n \n";
	else
	{
		std::cout<<" \nInput error : number of points in xhist dont match Jhist \n";
		return;
	}

    MatD poll_points,poll_to_eval;
    VecD xmin(dim);
    double Jmin;
    find_xmin(x,J,xmin,Jmin);
    poll_mads(xmin,dim,set0.spc,set0.xlow,set0.xupp,set0.nstage_poll,poll_points);
    int eval_size = 0;
    for(int i=0;i<dim+1;i++)
    {
    	if((newpoint(num_pts,dim,x,poll_points[i]))&&(bndchk(dim,poll_points[i],set0.xlow,set0.xupp)))
    	{
    		poll_to_eval.push_back(poll_points[i]);
    		eval_size++;
    	}
    }

	set0.npts_poll += eval_size ; 
	set0.nstage_poll += 1 ; 

    
	print_matrix_to_txt_file("xnew.dat",poll_to_eval,eval_size,dim,xlog_index);
}

void refine(int dim,VecI &spc)		// FOR UNSUCCESSFUL SEARCH AND UNSUCCESSFUL POLL
{
	for(int i=0;i<dim;i++)
		spc[i] = (spc[i]-1)*4 + 1;
}



void coarsen(int dim,VecI &spc)		// FOR EITHER A SUCCESSFUL SEARCH OR A SUCCESSFUL POLL
{
	double coarsen_fact = 1;

	std::cout<<"SPC[0] BEFORE : "<<spc[0]<<std::endl;
	for(int i=0;i<dim;i++)
		spc[i] = (spc[i]-1)/coarsen_fact + 1;
	std::cout<<"SPC[0] AFTER COARSEN : "<<spc[0]<<std::endl;
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

//	PARAMETER INDEX - TOGGLE FOR 0 - NORMAL; 1 -LOG SCALE
	VecI xlog_index(dim,0);

//	THETA INDEX VECTOR
	VecI th_index(dim);

//  SET BOUNDS FOR THETA
	VecD thlow(dim_th);
	VecD thupp(dim_th);
	VecD theta(dim_th);

	set_th_index(th_index);
	set_th_bounds(thlow,thupp);

	FILE *fp;
	fp = fopen("BestJHist.dat","a");


	settings set0 ;

	int coarsen_flag=0;

//=============== DEFINE BOX CONSTRAINTS ======================

	set_xlog(xlog_index);	
	read_phase("smf.inp",set0);

	std::cout<<"Past the read phase \n";

	bool dim_bds = ((set0.xupp.size()!=DIM_1)||(set0.xlow.size()!=DIM_1));

    if(dim_bds)
    {
        throw std::runtime_error("Dimensions of upper bounds/lower bounds are incorrect \n");
    }



//=============== SMF PARAMETERS ==============================

	double conv_tol = 0.000001;
	double x_del_sum = 0;
	VecD x_delta(dim);


//============== TERMINATION CONDITION ======================

    for(i=0;i<dim;i++)
    {
            x_delta[i] = (set0.xupp[i] - set0.xlow[i])/(set0.spc[i]-1);
            x_del_sum += x_delta[i];
    }


	if(x_del_sum <= dim*conv_tol)
	{
		std::cout<<"Insert post-processing stuff here \n";

		terminate_opt(0);
		return 0;
	}
	std::cout<<"X_DEL_SUM = "<<x_del_sum<<std::endl;


//=============== READ CURRENT STATE IN SMF LOOP =============

	if (set0.phase==0)  //  FIRST EXECUTION OF COMMAND
	{

		std::cout<<"LHS initiated \n";
		num_pts = N_DIV1;
		lhs(dim, num_pts, set0, 1, x, xlog_index)	;	
		
		print_matrix_to_txt_file("xnew.dat",x,num_pts,dim,xlog_index);
		
		//  INITIALIZE THETA
		for(i=0;i<dim_th;i++)
			theta[i]= sqrt(thupp[i]*thlow[i]);
		set0.theta = theta ; 
		set0.phase = 1;
		write_phase("smf.inp",set0);

		return 0;
	}

	else if (set0.phase==1)	// LHS COMPLETED  -- DO SEARCH
	{
		std::cout<<" \nLHS completed. Searching \n";
		int num_pts,num_pts_y;
		import_matrix_from_txt_file("xnew.dat",x,num_pts,dim,xlog_index);
		import_vector_from_txt_file("Jnew.dat",J,num_pts_y);
		
		if(num_pts!=num_pts_y)
			std::cout<<"Diff in xnew and Jnew = "<<num_pts-num_pts_y<<std::endl;

		print_matrix_to_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
		print_vector_to_txt_file("Jhist.dat",J,num_pts);

		bool new_pt = search(dim,dim_th,thlow,thupp,th_index,set0,xlog_index);

		if(new_pt == true)
		{	

			set0.phase = 2 ; 
		}
		else		
		{			

			std::cout<<"Pre-evaluated point found from search. Do poll \n";

			poll(dim,set0,xlog_index);
			set0.phase = 3;
		}

	}

	else if (set0.phase==2)	// SEARCH COMPLETED -- CHECK IF POLL IS NEEDED
	{
		bool search_flag = succ_search(dim,set0,xlog_index);
		std::cout<<"Checking if search worked \n";
		if(search_flag == true)
		{
			std::cout<<"Search worked \n";
			coarsen_flag=1;
			bool new_pt = search(dim,dim_th,thlow,thupp,th_index,set0,xlog_index);
			if(new_pt == true)
				set0.phase = 2 ; 

			else		// SEARCH GIVES A PRE-EVALUATED GRID POINT AS MINIMUM
			{
				std::cout<<"A pre-evaluated point obtained. Do poll \n";			
				// DO POLL
				poll(dim,set0,xlog_index);
				set0.phase = 3;
			}

		}
		else
		{
			std::cout<<"Search failed. Polling \n";
			poll(dim,set0,xlog_index);
			set0.phase = 3;
		}

	}
	else if (set0.phase==3)	// POLL COMPLETED -- CHECK IF POLL WORKED
	{
		std::cout<<"Checking if poll worked \n";
		bool poll_flag = succ_poll(dim,set0,xlog_index);

		if(poll_flag == true)
		{
			std::cout<<"Poll worked. Do search \n";
			coarsen_flag=1;

			bool new_pt = search(dim,dim_th,thlow,thupp,th_index,set0,xlog_index);
			if(new_pt == true)
				set0.phase = 2; 

			else		// SEARCH GIVES A PRE-EVALUATED GRID POINT AS MINIMUM
			{			// DO POLL
				std::cout<<"Search failed. Back to poll  \n";
				poll(dim,set0,xlog_index);
				set0.phase = 3; 
			}
		}
		else
		{
			std::cout<<"Poll failed. Refine and search\n";
			refine(dim,set0.spc);

			bool new_pt = search(dim,dim_th,thlow,thupp,th_index,set0,xlog_index);
			if(new_pt == true)
				set0.phase = 2; 

			else		// SEARCH GIVES A PRE-EVALUATED GRID POINT AS MINIMUM
			{			// DO POLL
				std::cout<<"Search failed. Do poll\n";
				poll(dim,set0,xlog_index);
				set0.phase = 3; 
			}
		}
	


	}
	else if (set0.phase == 4) // RESTART CASE
	{

		std::cout<<"Restarted SMF! \n";

		bool new_pt = search(dim,dim_th,thlow,thupp,th_index,set0,xlog_index);
		if(new_pt == true)
			set0.phase = 2;
		else		// SEARCH GIVES A PRE-EVALUATED GRID POINT AS MINIMUM
		{			// DO POLL
			poll(dim,set0,xlog_index);
			set0.phase = 3;
		}


	}


	else
	{
		std::cout<<"Error. The value of phase is : "<<set0.phase<<"\n";
		terminate_opt(set0.phase);
	}

	if(coarsen_flag==1)
		coarsen(dim,set0.spc);


	// ------- WRITE PHASE FILE 

	write_phase("smf.inp",set0);


	// VecD spc_d(dim);
	// for(int j=0;j<dim;j++)
	// 	spc_d[j] = (double)spc[j];

	// print_vector_to_txt_file("spc.dat",spc_d,dim);

	if(set0.phase>=2)
	{
		VecD xmin_final;
		double Jmin_final;
		int num_pts,num_pts_y;
		import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
		import_vector_from_txt_file("Jhist.dat",J,num_pts_y);

		find_xmin(x,J,xmin_final,Jmin_final);
		write_hist(fp,xmin_final,Jmin_final,xlog_index);

	    std::cout<<"Jmin till now : "<<Jmin_final<<" at xmin : \n ";
       	for(int p =0 ; p < dim ; p++)
           	std::cout<<xmin_final[p]<<"  ";
       	std::cout<<std::endl;
	}


	fclose(fp);
}

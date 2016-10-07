
#include "krig.h"

/*
------------------------------------------
	MADS BASED POLLING FUNCTION
------------------------------------------

INPUTS : x 	- Point to be polled
 	 n 	- Dimension of x
	 delta 	- Polling distance
	 spc	- Number of points along a cartesian grid direction
	 n_poll - (new) Counter for number of polls done till now, incremented at the end (initially = 0)
	 amin	- Minimum parameter value
	 amax 	- Maximum parameter value	 

OUTPUTS : poll_pts - Points generated around the given x 
	( Empty memory taken as input )


Polling routine re-adapted to match algorithm in Audet and Dennis' 2006 paper detailing the 
practical implementations of LTMADS for n+1 directions

*/

void poll_mads(VecD x, int n, VecI spc, VecD amin,VecD amax, int &n_poll, MatD &poll_pts)
{

	int r,i,j;
	int L[n-1][n-1],B[n][n];
	VecI tot(n,0);

	int D[n][n+1];

	VecI permr,permc; 
	VecD interv(n);
	// Generate seed here 

	std::srand(std::time(0));

	// Step 1 : Select a row number randomly  

	int base_row_id = randn(n);

	// Step 2 : Generate first base direction 
	// 2a  : Fill in (+/-) max direction size at the selected row number

	VecI base_row(n); 


	int max_dsize = pow(2,n_poll); 
	int dir = 2*(randn(2)) - 1; 

	int base_val = dir*max_dsize ; 




	// 2b  : Fill in complete row with {-2^npoll +1 , ... , 2^npoll -1 }


	for(i=0; i<n ;i++)
	{
		base_row[i] = randn( 2*(max_dsize-1) + 1 ) - max_dsize + 1;  
	}

	base_row[base_row_id] = base_val ; 


    /*DEBUG STEP
    //Print base row id
    std::cout<<"\nThe base row id is : "<<base_row_id<<" \n";
    // Print base_row vector
    std::cout<<"\nThe base row is :\n";
    for(i=0;i<n;i++)
    {
  		std::cout<<base_row[i]<<"\t";
    }
    std::cout<<std::endl;
    //END DEBUG*/

	// Step 3 : Generate (N-1) dim Lower triangular matrix

	for(i=0;i<n-1;i++)
	{
		dir = (randn(2))*2 - 1;
		L[i][i] = dir*max_dsize;

		for(j=0;j<i;j++)
		{
			L[i][j] = randn( 2*(max_dsize-1) + 1 ) - max_dsize + 1  ;
			L[j][i] = 0;
		}
	}


	// Step 4 : Generate n dimensional basis B 

	// 4a : Generate permutation vectors for row and columns of B

	// Row vector permutation lies in {0,..,n-1}\{base_row_id}. This will be used in the next step. 
	// Col vector permutation lies in {0,..,n}. This will be used right at the very end


	// This fills an ordered list of contents for both permutation vectors
	for(i =0;i<n;i++)
	{
		if(i!=base_row_id)
			permr.push_back(i);

		permc.push_back(i);

	}
	permc.push_back(n);

	// The ordered lists are shuffled here
	std::random_shuffle(permr.begin(),permr.end(),randn);
    std::random_shuffle(permc.begin(),permc.end(),randn);


    // 4b : Generate B as a combination of row permuted L and base row at location set by base_row_id

    // Set a vector iterator to start at the beginning of permr
    std::vector<int>::iterator v;

    v = permr.begin();

    for(i=0;i<n;i++)
    {
    	if(i<n-1)
    	{
    		for(j=0;j<n-1;j++)
    		{
    			B[*v][j] = L[i][j];
    		}
    		v++;

    		B[base_row_id][i] = 0 ;
    	}
    	B[i][n-1] = base_row[i];
    }


    /*DEBUG STEP
    // Print B matrix
    std::cout<<"\nThe basis matrix is :\n";
    for(i=0;i<n;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		std::cout<<B[i][j]<<"\t";
    	}
    	std::cout<<std::endl;
    }
    std::cout<<std::endl;
    //END DEBUG*/


    // Step 5: Generate +ve basis D 

    // Fill permuted columns of B into D
    v = permc.begin();

    for(j=0;j<n;j++)
    {
    	for (i = 0; i < n ; i++)
    	{
    		tot[i] = tot[i] + B[i][j];
    		D[i][*v] = B[i][j];
    	}
    	v++;
    }

    // Fill in the column = -sum(columns ofB) to complete positive basis
    for(i=0;i<n;i++)
    {
    	D[i][*v] = -tot[i];
    }

    //DEBUG STEP
    /* Print D matrix
    std::cout<<"\nThe positive basis matrix is :\n";
    for(i=0;i<n;i++)
    {
    	for(j=0;j<n+1;j++)
    	{
    		std::cout<<D[i][j]<<"\t";
    	}
    	std::cout<<std::endl;
    }
    //END DEBUG*/

    // Step 6: Convert D to grid units 

    // 6a : Compute minimum grid size for each dimension
    for(i=0;i<n;i++)
    {
    	interv[i] = (amax[i] - amin[i])/(spc[i]-1.0);
    }

    // 6b : Compute poll points
    for(j=0;j<n+1;j++)
	{	
		VecD temp(n,0);
		for(i=0;i<n;i++)
		{
			temp[i] = x[i] + D[i][j]*interv[i];
		}
		poll_pts.push_back(temp);
	}

}

//====================================================================


/* FOR STANDALONE TESTING- MAIN : 


int main()
{

	// Generate seed
//	std::srand(2);

	int n = 3;
	VecD x(n),amin(n),amax(n);
	double delta;
	int i,j;
	std::vector<int> spc(n);

	int stage = 3;

	MatD poll_pts;
//	poll_pts = (double **)std::malloc((n+1)*sizeof(double *));
//	for(i=0;i<n+1;i++) { poll_pts[i] = (double *)std::malloc(n*sizeof(double));}

//	poll_pts[0][0] = 1.0;


	for(i=0;i<n;i++)
        {
                x[i]=50;
		amin[i] = 0.01;
		amax[i] = 100;
		spc[i] = 101;
        }

	std::cout<<"Before poll \n";

	poll_mads(x,n,spc,amin,amax,stage,poll_pts);

	delta = (amax[0] - amin[0])/(spc[0]-1);
	std::cout <<"Delta : "<<delta<<" \n";

	std::cout<<"Original point : \n";
	for(i=0;i<n;i++){ std::cout<<x[i]<<"   "; }
	std::cout<<"\n Poll Points : \n ";

	for(i=0;i<n+1;i++)
	{
		std::cout<<" Point "<<i<<" \n "; 
		for(j=0;j<n;j++)
		{	std::cout<<poll_pts[i][j]<<" ";
		}
		std::cout<<"\n";
	}


}	


//END MAIN */

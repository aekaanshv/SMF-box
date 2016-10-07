#include "krig.h"

/*
------------------------------------------
	MADS BASED POLLING FUNCTION
------------------------------------------

INPUTS : x 	- Point to be polled
 	 n 	- Dimension of x
	 delta 	- Polling distance
	 spc	- Number of points along a cartesian grid direction
	 amin	- Minimum parameter value
	 amax 	- Maximum parameter value	 

OUTPUTS : poll_pts - Points generated around the given x 
	( Empty memory taken as input )

*/


//int randn(int n) { return std::rand()%n; }

void poll_mads(VecD x, int n, VecI spc, VecD amin,VecD amax, MatD &poll_pts)
{
//	double poll_pts[n+1][n];
	VecD d1(n), interv(n),delta(n);
	double D[n+1][n],Dperm[n][n],tot;
	int i,r,j;
	VecI permr;
	VecI permc;
	for(i=0;i<n;i++)
	{
		delta[i] = (amax[i]- amin[i])/(spc[i]-1.0);
		d1[i] = sqrt(1.00/delta[i]);
	}
// CONSTRUCT LOWER TRIANGULAR BASIS
	int lim=0;
	for(i=0;i<n;i++)
	{
		r = (randn(2))*2 - 1;
		D[i][i] = r*d1[i];
		lim = floor(d1[i]-1);
		if(lim<=0)
			lim = 0;
		for(j=0;j<i;j++)
		{
			D[i][j] = randn(2*lim+1) - lim ;
//			D[i][j] = randn(2*n+1) - n;
			D[j][i] = 0;
		}

// CONSTRUCT PERMR AND PERMC
		permr.push_back(i);
		permc.push_back(i);
	}

              
	std::random_shuffle(permr.begin(),permr.end(),randn);
        std::random_shuffle(permc.begin(),permc.end(),randn);

	std::vector<int>::iterator v;

	v = permr.begin();

// PERMUTE ROWS
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++){ Dperm[i][j] = D[*v][j]; }
		v++;
	}

	v = permc.begin();

// PERMUTE COLUMNS
        for(i=0;i<n;i++)
        {
                for(j=0;j<n;j++){ D[j][i] = Dperm[j][*v]; }
		v++;
        }

//      std::cout<<"things permuted \n";


// D HAS THE PERMUTED MATRIX. ADD THE (n+1)th ROW

	for(j=0;j<n;j++)
	{
		tot =0;
		for(i=0;i<n;i++) { tot = tot + D[i][j];}
		D[n][j] = - tot;

// COMPUTE INTERVAL
		interv[j] = (amax[j] - amin[j])/(spc[j]-1.0);
	}
/* DEBUG IF PERM MATRIX IS FORMED CORRECTLY
	for(i=0;i<n+1;i++)
	{
		for(j=0;j<n;j++)
			std::cout<<D[i][j]<<"  ";
		std::cout<<std::endl;
	}

*/ //END DEBUG


// COMPUTE POLL POINTS

//        std::cout<<"before poll point calc \n";


	for(i=0;i<n+1;i++)
	{	
		VecD temp(n,0);
		for(j=0;j<n;j++)
		{
			temp[j] = x[j] + D[i][j]*interv[j];
		}
		 poll_pts.push_back(temp);
	}

//        std::cout<<"function completed \n";
//	return 0;	
}

//====================================================================


/* FOR STANDALONE TESTING- MAIN : 

int main()
{
	int n = 5;
	VecD x(n),amin(n),amax(n);
	double delta;
	int i,j;
	std::vector<int> spc(n);

	MatD poll_pts;
//	poll_pts = (double **)std::malloc((n+1)*sizeof(double *));
//	for(i=0;i<n+1;i++) { poll_pts[i] = (double *)std::malloc(n*sizeof(double));}

//	poll_pts[0][0] = 1.0;


	for(i=0;i<n;i++)
        {
                x[i]=50;
		amin[i] = 0.01;
		amax[i] = 100;
		spc[i] = 5;
        }

	std::cout<<"Before poll \n";

	poll_mads(x,n,spc,amin,amax,poll_pts);

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


END MAIN */

#include "krig.h"
void push_to_grid(VecD &x, VecD xmin, VecD xmax, VecI spc)
{
	int n = x.size();
	int floor;
	double delta;
	for (int i=0;i<n;i++)
	{
		delta = (xmax[i]-xmin[i])/(spc[i]-1);
		floor = (int)(x[i]/delta);
		if (x[i]-floor*delta <= delta/2)
			x[i] = floor*delta;
		else
			x[i] = (floor+1)*delta;
	}	

}

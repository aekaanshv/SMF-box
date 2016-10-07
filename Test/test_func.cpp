#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <ctime>
#include <algorithm>

#define VecD std::vector<double>

double blackbox(VecD x)

{
// ENTER BLACK BOX FUNCTION HERE
        int n = x.size();
        double tot = 0;
	double tot2 = 0;
	double shift = 5;

 //  SPHERE 

        for(int i=0;i<n;i++)
                tot += ((x[i]-shift)*(x[i]-shift) ) ; 

        return tot;
}

void import_vector_from_txt_file(const char* filename_Y, VecD & v, int& num_pts)
{
    std::ifstream file_Y;    // Input file stream class to operate on files
    std::string line;        // String class to operate on strings
    file_Y.open(filename_Y);
    if (file_Y.is_open())
    {
        num_pts = 0;
        while(std::getline(file_Y, line))
        {
                std::istringstream is( line );
                // is = Input stream class to operate on strings
                double n;
                while( is >> n )
                {
                        v.push_back( n );
                        num_pts++;
                }
        }
//      std::getline(file_Y, line);
//      num_pts = ReadNumbers(line,v);
        file_Y.close();
    }
    else
    {
        std::cout << "file open failed";
    }

}



int main()
{
	int dim;
	VecD inp;
	import_vector_from_txt_file("c.txt",inp,dim);

	double value = blackbox(inp);

	FILE *fp;

	fp = fopen("sse.txt","w");
	fprintf(fp,"%f\n",value);
	fclose(fp);	
}

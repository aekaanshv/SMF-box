/*         
            MATRIX PARSER
Adapted from : http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/

*/

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

// USAGE :     import_matrix_from_txt_file("x.txt",v,rows,cols);


/* DEBUG : FOR STANDALONE USAGE

int ReadNumbers( const string & s, VecD v );
void import_matrix_from_txt_file(const char* filename_X, MatD & v, int& rows, int& cols);

int main()
{
    vector<vector <double> > v;
    int rows=0;
    int cols=0;

    import_matrix_from_txt_file("x.txt",v,rows,cols);
}
*/ // END DEBUG


int ReadNumbers( const string & s, VecD first ) 
{
    istringstream is( s );  // Input stream class to operate on strings
    double n;
    while( is >> n ) 
    {  // The >> is an extraction operator and parses input and places into the specified variable.
                        // It has a boolean return type
        first.push_back( n ); // Adds new element to the end of the vector.
    }
    return first.size();    // Returns the size of the final vector.
}



void import_matrix_from_txt_file(const char* filename_X, MatD v, int & rows, int & cols)
{
    vector <double> first_line;  // Vector that is read in and added into first line 
    ifstream file_X;    // Input file stream class to operate on files
    string line;        // String class to operate on strings
    int MAX_INT = 32766;
    file_X.open(filename_X);
    if (file_X.is_open())
    {
        int i=0;
        getline(file_X, line);  // Get line from the file
        
        
        cols =ReadNumbers( line, first_line ); // Read first line 
        
	    v.push_back(first_line);
        
        for ( i=1;i<MAX_INT+1;i++)
        {         // Read lines continuously till we hit the end of the input or MAX_INT
            if ( getline(file_X, line) == 0 ) break;
	        VecD temp;  	// Temporary vector to be read 
            ReadNumbers( line, temp );     // Read numbers from the line
            v.push_back(temp);
        }
        
        rows=i; 
        if(rows >MAX_INT) cout<< "Number of rows must be smaller than MAX_INT";
        
        file_X.close();
    }
    else
    {
        cout << "file open failed";
    }

}


void print_matrix_to_txt_file(const char* filename_X, MatD v, int num_pts, int dim)
{
    FILE *file_X;    // Input file stream class to operate on files
    file_X.fopen(filename_X);
    for(int i=0;i<num_pts;i++)
    {    
        for(int j=0;j<dim;j++)
            fprintf(file_X,"%f  ",v[i][j]);
        fprintf(file_X,"\n");
    }
        
    fclose(file_X);

}


void print_vector_to_txt_file(const char* filename_Y, VecD v, int dim)
{
    FILE *file_X;    
    file_X.fopen(filename_X,"w");

    for(int j=0;j<dim;j++)
        fprintf(file_X,"%f  ",v[j]);
        
    fclose(file_X);
}
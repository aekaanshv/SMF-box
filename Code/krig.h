// INCLUDE ALL HEADERS HERE. 
// DEFINE THE MOST COMMONLY USED COMMANDS - SHORTEN NAMESPACES 

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <ctime>
#include <algorithm>
#include <cctype>
#include <functional> 
#include <locale>

//using namespace std;

// DEFINE STATEMENTS
#define VecI std::vector<int>
#define VecD std::vector<double>
#define MatD std::vector< std::vector<double> >
#define MatI std::vector< std::vector<int> >

//DEFINE MODEL STRUCTURE

struct model
{
public:
    int num_pts,dim, dim_th;
    MatD x0,L;
	VecI th_index;
    VecD y0, theta;
};

//DEFINE SETTINGS STRUCTURE
struct settings
{
public:
	// if the value of phase is :
	// 0 --- START
	// 1 --- LHS 
	// 2 --- SEARCH 
	// 3 --- POLL
	// 4 --- RESTART (SEARCH/POLL STATS NOT USED IN THIS SETTING)
	int phase = -1; 

	VecI spc;

	VecD xupp, xlow;
	VecD theta;
	char* hist_id = NULL ; 

	int npts_lhs = 0 ;
	int npts_search = 0 ; 
	int npts_poll = 0 ; 

	int nstage_search = 0 ;
	int nstage_poll = 0 ; 

	int succ_search = 0 ; 
	int succ_poll = 0 ;
};

// FUNCTION PROTOTYPES - RANDOM NUMBER

int randn(int n); 

double randd (double lower_bound, double upper_bound);


// FUNCTION PROTOTYPES - FILE I/O

void read_phase(const char* filename_p, settings &);

void write_phase(const char* filename_p, settings );

void import_vector_from_txt_file(const char* filename_Y, VecD & v, int& num_pts);

void import_matrix_from_txt_file(const char* filename_X, MatD & v, int& rows, int cols, VecI xlog_index);

void print_matrix_to_txt_file(const char* filename_X, MatD v, int num_pts, int dim, VecI xlog_index);

void print_vector_to_txt_file(const char* filename_Y, VecD v, int dim);

void write_hist(FILE *fp, VecD x, double J, VecI xindex);

// FUNCTION PROTOTYPES - LHS

void lhs(int , int , settings &, bool , MatD &,VecI ) ;

// FUNCTION PROTOTYPES - KRIGING

bool Cholesky(int d,MatD R,MatD &L);

void expcov(int num_pts,int dim,int dim_th,MatD x,VecD x0,VecI th_index,VecD theta,VecD &covx);

void chol_sol(int ,MatD , VecD , VecD &);

double mean_cal(int num_pts, MatD L, VecD y0);

double sigma2_cal(int num_pts, MatD L, VecD y0, double mean);

double cal_det_tri(MatD L,int n);

double loglikelihood(model,double const *);

void predictor(model, int, MatD , VecD &);

double msecheck(model );

void push_to_grid(VecD &, VecD , VecD , std::vector<int> );

// FUNCTION PROTOTYPES - POLLING

void poll_mads(VecD , int , std::vector<int> , VecD ,VecD , int &, MatD &);

// FUNCTION PROTOTYPES - MODEL

void model2param(model,int &, int&, int &, MatD &, VecD &,VecI &,VecD &, MatD &);

bool model2param_wchol(model,int &, int &, int &, MatD &, VecD &,VecI &, VecD &, MatD &);

bool param2model(model &,int,int,int,MatD,VecD,VecI,VecD);

// FUNCTION PROTOTYPES - OPTIMIZATION

//void genpattsearch(int dim,VecD amin,VecD amax, std::vector<int> spc, VecD &thmax, double &Jmax, model S, double (*feval)(model));

double pred_fn(model S, double const *x);

void nelmin (model,double fn(model,VecD),  int, VecD &, VecD &, double *, double, VecD , int , int, int *, int *, int * );

double * cm_optimize(int dim, model S, double(*pFun)(model, double const *), VecD xlow, VecD xupp, int number_of_restarts, double increment_factor_for_population_size, VecD inguess, VecD in_stddev, double tol);

bool bndchk(int dim, VecD points, VecD amin, VecD amax);


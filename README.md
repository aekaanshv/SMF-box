# SMF-box
C++ code for Box constrained SMF 

Directories: 

* Code- Code for the SMF 
* Test- A test function, input file and example coupling script for the SMF

Settings that the user may want to modify from within the code:

* Dimension of the problem ( main.cpp, macro name DIM_1 )
* Number of LHS points ( main.cpp, macro name N_DIV1 )
* Log scale parameters ( main.cpp, variable name xlog_index. Follow instructions in code )
* Kriging regularization parameter ( likelihood.cpp, variable name k_penalty )
* Number of points tested per Hyper-parameter optimization  (main.cpp, variable name n_test_pts )

Compilation: 

* Downlaod repository and go into the code directory
* Type make to compile the code



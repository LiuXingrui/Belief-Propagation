#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>

#ifndef BP
#define BP
#include"BP.h"
#endif


/** @brief Belief-propagation functions
   
    @author Xingrui Liu <xliu@ucr.edu>
    @date August 2022
    */

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


/*! \relates nodes
 *Perform the matrix reduction for H and D, for stabilizer codes, if H=Hx, D=Hz, I call it D because it represents the degeneracy

\param A stores the columns of D which are added.
\param num_of_row_reduction is the number of how many rows are reduced.
 */
//
void Mat_Trans(GF2mat& H, GF2mat& D, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug=0);


void row_reduction(int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<vector<int>>& A, int debug=0);

void Add_cols (GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug=0);

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
void Mat_Trans(const GF2mat& H, const GF2mat& D, vector<double> &K, vector<double>& K_tilde, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug=0);


void row_reduction(const int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<double> &input_K,vector<double>&output_K,vector<vector<int>>& A, int debug=0);

void construct_B_tau(vector<double> &B_tau,int w, vector<double>& input_K, vector<int> Ai,vector<vector<int>>&  Tau);

void Add_cols(GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug=0);

void K_trans(vector<double>& input_K, vector<double> &output_K, vector<int>& Ai, int w,vector<vector<int>>& b);

void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug,vector<vector<int>>& b);

void construct_BF(GF2mat& D,vector<int> &Ai, int w,int c, GF2mat& DB, GF2mat& F,vector<vector<int>>& b,int n,int debug=0);

GF2mat merge_mat_vert(const GF2mat &up,const GF2mat &bottom, int debug=0);

bool new_decoder1(GF2mat& H2,GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,double& num_iter, int lmax,int& max_fail,int&syn_fail,int debug, vec &LR,int rankH,int options);

bool new_decoder2(GF2mat& H,GF2mat& G,GF2mat &H_tilde, GF2mat &H_tilde_star,const vector<vector<int>>& A,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,const int&  num_row_red,double& num_iter, int lmax,int& max_fail,int&syn_fail,int& trans_suc,int & no_error,int debug, vec &LR,int rankH,int options);
void inverse_trans(GF2mat& output_e,const GF2mat& output_e_tilde,const vector<vector<int>>& A);

int Weight(const bvec &cw);

double energy (bool empty_G,const GF2mat& G, const GF2mat& eT, const vector<double> &K);
double ML_suc_rate (vec p,const GF2mat& D,const GF2mat& H, const GF2mat& H_tilde_star, vector<double>&K,int d,int num_large_wt_error);

void ML_decoder_verify(vec p,const GF2mat& H, const GF2mat D,const GF2mat& H_star, const GF2mat D_star,const GF2mat H_tilde, const GF2mat H_tilde_star,vector<double>&K,int max_num_cws,double& after_trans_suc_rate,double&suc_rate,double &after_trans_theoric_suc_rate,double& min_wt_suc_count,int debug=0);

void ML_decoder(const GF2mat& input_e, const GF2mat& D,const GF2mat& L,vector<double>& K,GF2mat& output_e);

void algebraic_decoder(const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e);

void get_logical(const GF2mat& H, const GF2mat& D,const GF2mat& H_star, GF2mat& L);
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
void Mat_Trans(const GF2mat& H, const GF2mat& D, vector<double> &K, vector<double>& K_tilde, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug=0);


void row_reduction(const int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<double> &input_K,vector<double>&output_K,vector<vector<int>>& A, int debug=0);

void construct_B_tau(vector<double> &B_tau,int w, vector<double>& input_K, vector<int> Ai,vector<vector<int>>&  Tau);

void Add_cols(GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug=0);

void K_trans(vector<double>& input_K, vector<double> &output_K, vector<int>& Ai, int w,vector<vector<int>>& b);

void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug,vector<vector<int>>& b);

void construct_BF(GF2mat& D,vector<int> &Ai, int w,int c, GF2mat& DB, GF2mat& F,vector<vector<int>>& b,int n,int debug=0);

GF2mat merge_mat_vert(const GF2mat &up,const GF2mat &bottom, int debug=0);

bool new_decoder1(GF2mat& H2,GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,double& num_iter, int lmax,int& max_fail,int&syn_fail,int debug, vec &LR,int rankH,int options);

bool new_decoder2(GF2mat& H,GF2mat& G,GF2mat &H_tilde, GF2mat &H_tilde_star,const vector<vector<int>>& A,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,const int&  num_row_red,double& num_iter, int lmax,int& max_fail,int&syn_fail,int& trans_suc,int & no_error,int debug, vec &LR,int rankH,int options);
void inverse_trans(GF2mat& output_e,const GF2mat& output_e_tilde,const vector<vector<int>>& A);

int Weight(const bvec &cw);

double energy (bool empty_G,const GF2mat& G, const GF2mat& eT, const vector<double> &K);
double ML_suc_rate (vec p,const GF2mat& D,const GF2mat& H, const GF2mat& H_tilde_star, vector<double>&K,int d,int num_large_wt_error);

void ML_decoder_verify(vec p,const GF2mat& H, const GF2mat D,const GF2mat& H_star, const GF2mat D_star,const GF2mat H_tilde, const GF2mat H_tilde_star,vector<double>&K,int max_num_cws,double& after_trans_suc_rate,double&suc_rate,double &after_trans_theoric_suc_rate,double& min_wt_suc_count,int debug=0);

void ML_decoder(const GF2mat& input_e, const GF2mat& D,const GF2mat& L,vector<double>& K,GF2mat& output_e);

void algebraic_decoder(const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e);

void get_logical(const GF2mat& H, const GF2mat& D,const GF2mat& H_star, GF2mat& L);

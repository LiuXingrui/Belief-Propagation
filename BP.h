#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>


using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

class  nodes{
public:
  ivec neighbors;
  int degree;


};

//find the neighbors of check nodes 
void initialize_checks (const GF2mat &H, nodes checks[], int & E);

//find the neighbors of variable(error) nodes
void initialize_errors(const GF2mat &H, nodes errors[]);

//initialize massages to all-one matrices
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H_sparse);

//parallel schedule for classical BP
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);

//serial schedile
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);
  
//update vi to cj massage:
inline void update_ci_to_vj(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int i,int j,bin s);
inline void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr,double alpha=1);


//clasical decoding
int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,vec& pv,int debug);

// distance between 2 codewords
int distance(const GF2mat &output_e, const GF2mat &real_e,int n);

//read a parity check matrix
GF2mat read_matrix (int& n,int &r, string & file_name);

//write a matrix to a file
void write_matrix(string file_name, GF2mat &H);

//merge 2 matrices horizontally 
GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right);

  
//the bsc error channel
int error_channel(GF2mat &cw, const vec &p);

// error channel for fixed weight errors
void error_channel2(GF2mat &error, int wt);

//weight of a column vector
int s_weight(const GF2mat &s);

//weight of a row vector
int weight(const GF2mat &cw);

//get the error rate distribution
void pro_dist(double pmin,double pmax, vec& pv);

// for quantum decoding
bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,int& max_fail, int&syn_fail,  int debug, vec &LR,int rankH,int &OSD_suc,double alpha,double lambda);

  //ordered statistical decoder
void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat& output_e);

  	
//get a column permutation matrix
GF2mat col_permutation_matrix(ivec & perm);
GF2mat col_permutation_matrix_s(ivec & perm);

//quantum parallel schedule
void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha);

//quantum serial schedule
void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha);

//print the positions of the errors for checkerboard codes.
void err_pos2(const GF2mat &error);

// the rank of GF2mat matrix
int GF2mat_rank(const GF2mat& H_s);

//convert dense GF2mat to sparse
void dense_to_sparse(GF2mat &G,GF2mat& G_s);

void err_pos1(const GF2mat &error);

//get the generator matrix
GF2mat get_gen(const GF2mat &H);

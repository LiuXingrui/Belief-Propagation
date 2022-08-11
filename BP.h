#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>


/** @brief Belief-propagation functions

   
    @author Xingrui Liu <xliu@ucr.edu>
    @date August 2022
    */

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

class  nodes{

public:
  ivec neighbors; /**< store the neighbors of this node */
  int degree; /**< degree of this node */
};

/*! \relates nodes
 * Initilize the check nodes, for each check, get its neighbors and degree. 

\param E stores the number of edges in this graph.
\param H is the parity check matrix.
 */
void initialize_checks (const GF2mat &H, nodes checks[], int & E);


/*! \relates nodes
 *Initilize the variable nodes.
 */
void initialize_errors(const GF2mat &H, nodes errors[]); 

/*! \relates nodes
 *Initilize the messages , set   mcv(i,j) and   mvc(i,j) to 1 if H(i,j)=1. 

\param mcv stores the messages from check nodes to variable nodes. 
\param mvc stores the messages fom variable nodes to check nodes.
 */
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H);

/*! \relates nodes

 *Do a parallel iteration for classical BP.
\param c is the number of check nodes.
\param v is the number of varaible nodes.
\param output_e is the result error vector.
 */
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e); 

/*! \relates nodes

 *Do a serial iteration for classical BP.
\param c is the number of check nodes.
\param v is the number of varaible nodes.
\param output_e is the result error vector.
 */
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);
  
/*! \relates nodes
 *  Update the message c_i to v_j, that is mcv(i,j).

\param s is syndrome(i).
 */
inline void update_ci_to_vj(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int i,int j,bin s);

/*! \relates nodes
 * Update the message v_j to c_i, that is mvc(i,j).

\param ipr =(1-p)/p.
\param alpha is for larger step size.
 */
inline void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr,double alpha=1);


/*! \relates nodes
 *Generate a error vector, decode it by classical BP, output 1 if secceeds, output 0 if fails.

 *\param v is the number of variable nodes.
 *\param c is the number of check nodes.
\param H is the parity check matrix.
 *\param num_iter stores the total iterations for successful decoding. num_iter/num_of_suc_decoding=averge iterations for seccessful decoding.
 * \param lmax is the maximum number of iterations.
 * \param er stores the number of bit errors after decoding, er/num_of_cws= bit error rate afer decoding.
 *\param debug debug=1 for parallel schedule, debug=0 for serial schedule.
 */
int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,double p,int debug); 

/*! \relates nodes
 * Output the distance between output_e and real_e.

\param n is the length of the vectors.
 */
int distance(const GF2mat &output_e, const GF2mat &real_e,int n);

/*! \relates nodes
 *  Read a parity check matrix from file_name, return it

\param  n stores the columns of the martrix
\param  r stores the rows of the matrix.
 */
GF2mat read_matrix (int& n,int &r, string & file_name);

/*! \relates nodes
 *Write a parity check matrix H in file_name.
 */
void write_matrix(string file_name, GF2mat &H);

/*! \relates nodes
 * Merge 2 matrices horizentolly, return a matrix=(left,right).
 */
GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right);

/*! \relates nodes
 * Classical BSC channel with error rate p, input codeword cw, return  weight of the error.
 */
int cla_error_channel(GF2mat &cw, double p);

/*! \relates nodes
 *Quantum BSC channel, e(i)=1 with probalility p(i). Return  weight of the error.
 */
int error_channel(GF2mat &cw, const vec &p);

/*! \relates nodes
 *A quantum error channel which generates errors with weight wt.
 */
void error_channel2(GF2mat &error, int wt);

/*! \relates nodes
 *Return the weight of a column vector.
 */
int s_weight(const GF2mat &s);

/*! \relates nodes
 *Return the weight of a row vector.
 */
int weight(const GF2mat &cw);

/*! \relates nodes
 *Get the probability distribution pv, pv(i) is randomely choosen betweem pmin and pmax.
 */
void pro_dist(double pmin,double pmax, vec& pv);

/*! \relates nodes
 *Generate a error vector, decode it by quantum BP, return true if secceeds, return false if fails.

 *\param H is the parity check matrix.
\param G is the matrix for corresponding logical operators.
 *\param pv is the error rate distribution for error channel.
\param pv_dec is the error rate distribution for decoding.
 *\param num_iter is for calculating average number of iterations.
 *\param lmax is maximum number of iterations.
 *\param wt if wt=0, use error_channel, if wt>=1, use error_channel2.
 *\param max_fail stores the number of failures such that reaches the maximum iterations.
 \param syn_fail stores the number of failures such that output_e+real_e=a logical operator.
 \param debug is for bitwise options, see readme.md
 \param LR is the likelyhood-ratios.
 *\param rankH is the rank of H.
 *\param OSD_suc stores the number of successful decoding by OSD.
 *\param alpha is for larger step size.
 * \param lambda has not been used yet, just set it be 1 is ok.
 */
bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,int& max_fail, int&syn_fail,  int debug, vec &LR,int rankH,int &OSD_suc,double alpha,double lambda);

/*! \relates nodes
 *Ordered statistical decoder, input LR, get output_e.
 */
void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat& output_e);

  	
/*! \relates nodes
 *Generate a column permutation matrix by permutation vector perm.
 */
GF2mat col_permutation_matrix(ivec & perm);
GF2mat col_permutation_matrix_s(ivec & perm);

/*! \relates nodes
 *Do a quantum parallel iteration.
 */
void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha);

/*! \relates nodes
 *Do a quantum serial iteration.
 */
void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha);

/*! \relates nodes
 *Show the patterns of error for checkerboard codes.
 */
void err_pos2(const GF2mat &error);

/*! \relates nodes
 * Return the rank of a GF2mat matrix
 */
int GF2mat_rank(const GF2mat& H_s);

/*! \relates nodes
 *Convert dense GF2mat matrix to sparse one.
 */
void dense_to_sparse(GF2mat &G,GF2mat_sparse& G_s);

/*! \relates nodes
 *cout i if error(i)==1.
 */
void err_pos1(const GF2mat &error);

/*! \relates nodes
 *Return the logical operators matrix for H.
 */
GF2mat get_gen(const GF2mat &H);

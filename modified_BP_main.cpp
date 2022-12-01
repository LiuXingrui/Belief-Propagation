#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <chrono>


#ifndef BP
#define BP
#include"BP.h"
#endif

#ifndef modified_BP
#define modified_BP
#include"modified_BP.h"
#endif


using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


int main(int argc, char **argv){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//For timing.
  GlobalRNG_randomize ();
  double p;
  long long int num_of_cws=0;
  long long int max_num_of_cws=1e18;//when use number_of_decoding_failure, set max_num_of_cws=a very large number
  string file_name;
  string file_name2;
  string data_file;
  int lmax;
  int wt;
  int channel;

  int d=-1;
  int options=0;

  int num_of_failed_cws=0;
  int max_failed_cws;
  int num_row_red; //number of row reduction
  int debug=0;

  if (argc!=11){cout<<" need 11 parameters: ./MBP  Hx_file Hz_file p  max_failed_cws/max_num_of_cws lmax data_file options channel num_of_row_reduction debug"<<endl;return 1;}
  //get the parameters: 
   
  file_name=argv[1];
  file_name2=argv[2];
  data_file=argv[6];
 
  istringstream argv3( argv[3] );	
  if ( argv3 >> p){}
  else
    {
      cout<<"pavg should be a double"<<endl;
      return 1;
    }
  
  istringstream argv4( argv[4] );

 	     
  if ( argv4 >> max_failed_cws){}
  else
    {
      cout<<"max_failed_cws should be an int"<<endl;
      return 1;
    }

  //cout<<max_failed_cws<<endl;
  if (max_failed_cws<0)
    {
      max_num_of_cws=-max_failed_cws;
      max_failed_cws=100000000;// when use max_num_of_cws, set max_failed_cws to a large number.
    }
 

  istringstream argv5( argv[5] );
  if ( argv5 >> lmax){}
  else
    {
      cout<<"lmax should be an int"<<endl;
      return 1;
    }

  
  istringstream argv7( argv[7] );
  if ( argv7 >> options){}
  else
    {
      cout<<"options should be an int"<<endl;
      return 1;
    }

  istringstream argv8( argv[8] );
  if ( argv8 >> channel){}
  else
    {
      cout<<"channel should be an int"<<endl;
      return 1;
    }

  istringstream argv9( argv[9] );
  if ( argv9 >>  num_row_red ){}
  else
    {
      cout<<"  num_of_row_reduction should be an int"<<endl;
      return 1;
    }


  istringstream argv10( argv[10] );
  if ( argv10 >>debug ){}
  else
    {
      cout<<"  debug should be 0 or 1"<<endl;
      return 1;
    }

  double num_iter=0.0; //for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
 
  int max_fail=0;//number of fails that reach maximum iterations
  int syn_fail=0;//number of fails that get the right syndrome but fails
  bool Hx_suc=false;

  int n1,n2,n,k,r1,r2,r;
  
  
  GF2mat Hx=read_matrix( n1,r1, file_name);
  GF2mat Hz=read_matrix( n2,r2, file_name2);
  //Hx=Hx.get_submatrix(0,0,r1-1,n1-1);
  //Hz=Hz.get_submatrix(0,0,r2-1,n2-1);
  //are the parity check matrices right?
  if (n1!=n2){cout<<"nx!=nz, the two matrices donot match"<<endl;return 1;}  
  n=n1;

  r=r1+r2;
  int rankx=GF2mat_rank(Hx);
  int rankz=GF2mat_rank(Hz);
  int rankH=rankx;
  
  // if (num_row_red>=rankz){cout<<"error:   num_of_row_reduction>=rankz";return 1;}
  
  k=n-rankx-rankz;
  GF2mat zero_mat1(r1,r2);
  
  if((Hx*Hz.transpose())==zero_mat1) {}
  else{cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl;return 1;}
  if (debug==1){cout<<"input H(Hx) is \n"<<Hx<<"\n input D(Hz) is\n"<<Hz<<endl;}
  GF2mat G=get_gen(Hx);
  //cout<<Gx*(Hx.transpose())<<endl;
  
  GF2mat H=Hx;
  GF2mat D=Hz;
  
  GF2mat H_tilde,D_tilde;
  vector<vector<int>> A;
  vector<double> K,K_tilde;
  
  int K_initial=0.5*log((1-p)/p);
  
  for (int i=0;i<n;i++) {K.push_back(K_initial);}

  Mat_Trans(H,D,K,  K_tilde, H_tilde, D_tilde, A, num_row_red,debug);
  
  GF2mat H_tilde_star=get_gen(H_tilde);  

  nodes  checks[H_tilde.rows()];//checks for Hx and Z errors
  nodes errors[H_tilde.cols()];


  int E1=0; //number of edges in Hx factor graph, but this parameter is not used in this prog

  //find the neighbourhoods of all nodes:


  initialize_checks (H_tilde, checks, E1);
  initialize_errors(H_tilde, errors);

  vec pv(n);
  vec LR(K_tilde.size());
  vec pv_dec(K_tilde.size());
  // int  num_of_suc_dec=0;
  for (int i=0;i<n;i++){pv(i)=p;}
  for (int j=0;j<K_tilde.size();j++){pv_dec(j)=1/(1+exp(2*K_tilde[j]));}    

  
  while (num_of_failed_cws<max_failed_cws&&num_of_cws<max_num_of_cws)
    {
      num_of_cws++;
      Hx_suc=new_decoder2(H, G,H_tilde, H_tilde_star, A, checks,errors,pv,pv_dec,  num_row_red, num_iter,  lmax, max_fail,syn_fail,debug, LR, rankH, options);
 
      // cout<<num_iter<<endl;
      if (Hx_suc==true)
	{       	 
	  num_of_suc_dec++;
	}
      else
	{
	  num_of_failed_cws++;
	}  
    }
   



  cout<<"p="<<p<<", there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a [["<<n<<", "<<k<<"]] code (decode x errors only)"<<endl;
   cout<<"average iterations:"<<endl;
 
   cout<<num_iter/num_of_suc_dec<<"\n\n"<<endl;
   cout<<"syn_fail="<<syn_fail<<endl;
   cout<<"max_fail="<<max_fail<<endl;
  
   // cout<<"num of zero errors is about "<<pow(p,n)*num_of_cws<<endl;
      
   ofstream myfile;
   myfile.open (data_file,ios::app);
 
   myfile << n<<" "<<d<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<p<<" "<<1.0*num_iter/num_of_suc_dec<<"  "<<num_of_suc_dec<<" "<<num_of_cws<<"  "<<syn_fail<<" "<<max_fail<<" "<<1.0*syn_fail/num_of_cws<<" "<<1.0*max_fail/num_of_cws<<endl;
  

   myfile.close();
   std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
   std::cout << "\n Run-time " << time_span.count() << " seconds.\n";

  
  return 0;
}

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
  int num_large_wt_error;
  string file_name;
  string file_name2;
  int debug;
  
  if (argc!=6){cout<<" need 5 parameters: ./MLSR  Hx_file Hz_file p  num_large_wt_error debug"<<endl;return 1;}
  //get the parameters: 
   
  file_name=argv[1];
  file_name2=argv[2];
  
 
  istringstream argv3( argv[3] );	
  if ( argv3 >> p){}
  else
    {
      cout<<"pavg should be a double"<<endl;
      return 1;
    }



  istringstream argv4( argv[4] );	
  if ( argv4 >>  num_large_wt_error){}
  else
    {
      cout<<" num_large_wt_error should be an int"<<endl;
      return 1;
    }

  istringstream argv5( argv[5] );	
  if ( argv5 >>  debug){}
  else
    {
      cout<<" debug should be an int"<<endl;
      return 1;
    }
 
 
  int n1,n2,n,k,r1,r2,r;
  
  
  GF2mat Hx=read_matrix( n1,r1, file_name);
  GF2mat Hz=read_matrix( n2,r2, file_name2);

  
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

  
  GF2mat H=Hx;
  GF2mat D=Hz;
  
  GF2mat H_tilde,D_tilde;
  vector<vector<int>> A;
  vector<double> K,K_tilde;
  
  double K_initial=0.5*log((1-p)/p);
  // cout<<"K_initial is "<<K_initial<<endl;
  
  for (int i=0;i<n;i++) {K.push_back(K_initial);}

  Mat_Trans(H,D,K,  K_tilde, H_tilde, D_tilde, A, D.rows(),debug);
  /*
  cout<<"at the end K_tilde:"<<endl;
   for (int ii=0;ii<K_tilde.size();ii++)
    {
      cout<<K_tilde[ii]<<" ";
    }
   cout<<endl;
  */
  GF2mat H_tilde_star=get_gen(H_tilde);
  if (H_tilde_star.rows()!=1){cout<<"error: H_tilde_star.rows()="<<H_tilde_star.rows()<<endl; return 1;}

  int d=sqrt(n);
  

  vec pv(n);
  for (int i=0;i<n;i++){pv(i)=p;}


  double suc_rate= ML_suc_rate (pv,D,H,  H_tilde_star,K,d, num_large_wt_error);

  //  cout<<4<<endl;

  cout<<"for [n="<<n<<", k=1, d="<<d<<"] surface code, \n raw bit suc rate (1-p)="<<1-p<<"\n estimated ML suc rate="<<suc_rate<<" \n using "<<num_large_wt_error<<" errors with wt>=d/2"<<endl;
   

  


   std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
   std::cout << "\n Run-time " << time_span.count() << " seconds.\n";

  
  return 0;
}

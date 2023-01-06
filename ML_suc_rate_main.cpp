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
  int num_cws;
  string file_name;
  string file_name2;
  string data_file;
  int debug;
  
  if (argc!=7){cout<<" need 6 parameters: ./MLSR  Hx_file Hz_file p  num_cws debug data_file"<<endl;return 1;}
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
  if ( argv4 >>  num_cws){}
  else
    {
      cout<<" num_cws should be an int"<<endl;
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
  if (H_tilde_star.rows()!=1){cout<<"H_tilde_star has "<<H_tilde_star.rows()<<"rows left, that is too much for calculating ML suc rate,so quit the prog"<<endl;return 0;}

  GF2mat H_star=get_gen(H);
  GF2mat D_star=get_gen(D);
  
 
  
  int d=sqrt(n);
  

  vec pv(n);
  for (int i=0;i<n;i++){pv(i)=p;}

  double after_trans_suc_rate,suc_rate,after_trans_theoric_suc_rate,min_wt_rate;
  //cout<<1<<endl;
  
  ML_decoder_verify(pv,H, D, H_star,  D_star,H_tilde,H_tilde_star,K,num_cws,after_trans_suc_rate,suc_rate,after_trans_theoric_suc_rate,min_wt_rate,debug);

  //  cout<<2<<endl;
  // double suc_rate= ML_suc_rate (pv,D,H,  H_tilde_star,K,d, num_large_wt_error);


  //  cout<<4<<endl;

  cout<<"for [n="<<n<<", k=1, d="<<d<<"] surface code, \n raw bit suc rate (1-p)="<<1-p<<"\n  ML suc rate="<<suc_rate<<" \n after_trans_ML_suc_rate="<<after_trans_suc_rate<<"\n after_trans_theoric_suc_rate ="<<after_trans_theoric_suc_rate<<" \n min_wt_suc_rate="<<min_wt_rate<<" \n using "<<num_cws<<" codewords"<<endl;

    ofstream myfile;
    myfile.open (data_file,ios::app);
    //  myfile<<"# 1-p   ML_suc_rate  after_trans_ML_suc_rate   after_trans_theoric_suc_rate"<<endl;
    myfile<<p<<"  "<<suc_rate<<"     "<<after_trans_suc_rate<<"      "<<after_trans_theoric_suc_rate<<"   "<<min_wt_rate<<endl;

  


   std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
   std::cout << "\n Run-time " << time_span.count() << " seconds.\n";

  
  return 0;
}

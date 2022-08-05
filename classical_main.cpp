#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <chrono>

#ifndef BP
#define BP
#include"BP.h"
#endif


#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;
//compile: g++   `pkg-config --cflags itpp` -o my_prog BP.cpp  BP1.cpp `pkg-config --libs itpp`
// format : ./my_prog H_file data_file  pmin pmax number_of_cordwords max_iteration debug
//the format for the parity check file: the first row is n n-k,  then the sparse matrix. Notice
// the first element stored in the file is 1 rathar than 0.

int main(int argc, char **argv){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

  GlobalRNG_randomize ();
  double pmax;
  double pmin;
  int num_of_cws;
  string file_name;
  string data_file;
  int lmax;
  int debug;
  
  if (argc!=8){
    cout<<" # of parameters!=7:  ./my_prog H_file data_file  pmin pmax number_of_cordwords max_iteration debug"<<endl;
    return 1;
  }
  //get the parameters: 
    
  file_name=argv[1];
  data_file=argv[2];
 
  istringstream argv3( argv[3] );	
  if ( argv3 >> pmin){}
  else
    {
      cout<<"pmin should be a double"<<endl;
      return 1;
    }

  istringstream argv4( argv[4] );
  if ( argv4 >> pmax){}
  else
    {
      cout<<"pmax should be a double"<<endl;
      return 1;
    }
	   
	     
  istringstream argv5( argv[5] );
  if ( argv5 >> num_of_cws){}
  else
    {
      cout<<"num_of_cws should be an int"<<endl;
      return 1;
    }
	        
  istringstream argv6( argv[6] );
  if ( argv6 >> lmax){}
  else
    {
      cout<<"lmax should be an int"<<endl;
      return 1;
    }
   istringstream argv7( argv[7] );
  if ( argv7 >> debug){}
  else
    {
      cout<<"debug should be an int"<<endl;
      return 1;
    }
  if (debug==0) {cout<<"parallel schedule"<<endl;}
  else {cout<<"serial schedule"<<endl;}
  double num_iter=0.0; //for calculate average iterations
  int num_of_suc_dec=0;// number of successfully decoded results


  //read the parity check matrix:
  int n;
  int &v=n;
  int r;
  int &c=r;
 
  GF2mat H=read_matrix ( n,r,  file_name);
  //cout<<H<<endl;
  nodes  checks[c];
  nodes  errors[v];
  int E=0;
  //find the neighbourhoods of all nodes:
  initialize_checks (H, checks,  E);
  initialize_errors(H, errors);
  int k=n-GF2mat_rank(H);
  int er=0;  //er is the number of bits that are wrong after decoding

  vec pv(n);   
  pro_dist( pmin,pmax, pv);
  
  for (int s=0;s<num_of_cws;s++)
    {     
      num_of_suc_dec= num_of_suc_dec+cla_decode( v,c,H, checks, errors, num_iter,  lmax, er,pv,debug);
    }
  
  cout<<"for p in ("<<pmin<<", "<<pmax<<"), there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a ["<<v<<", "<<k<<"] code"<<endl;
  // cout<<"and there are "<< n_valid_cws<<" decoding results are codewords"<<endl;
  cout<<"bit error rate after decoding is ";
  cout<<er/(1.0*n*num_of_cws)<<endl;
  // cout<< "number of zero errors:    "<<num_zero_cws<<endl;
  cout<<"average iterations:"<<endl;
  cout<<num_iter/num_of_suc_dec<<endl;



  double midp=(pmax+pmin)/2;
      
  ofstream myfile;
  myfile.open (data_file,ios::app);
  myfile << n<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<midp<<" "<<num_iter/num_of_suc_dec<<"  "<<num_of_suc_dec<<"  "<<er/(1.0*n*num_of_cws)<<endl;
  myfile.close();
  std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
  std::cout << "\n Run-time " << time_span.count() << " seconds.\n";

  return 0;
}



    

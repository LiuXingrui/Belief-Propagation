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



using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

//GlobalRNG_reset (1);
int main(int argc, char **argv){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  GlobalRNG_randomize ();
  double pmax;
  double pmin;
  int num_of_cws=0;
  int max_num_of_cws=2147483647;
  string file_name;
  string file_name2;
  string data_file;
  int lmax;
  int wt;
  int channel;

  int d=-1;
  int debug=0;
  double pavg;
  double range;
  double alpha;
  double lambda;
  int num_of_failed_cws=0;
  int max_failed_cws;
  double decode_p,decode_prange,decode_pmin,decode_pmax;
  //can input pmin and pmax, or the weight of error_vectors and dec_method, which is same_p or diff_p: 
  //diff p decode: pmin=wt/n*0.5, pmax=wt/n*1.5,   same p decode: p=wt/n
  if (argc!=14){cout<<" need 13 parameters: ./qBPx  Hx_file Hz_file pavg/wt range  max_failed_cws lmax data_file debug channel alpha decode_p decode_prange"<<endl;return 1;}
  //get the parameters: 
   
  file_name=argv[1];
  file_name2=argv[2];
  data_file=argv[7];
 
  istringstream argv3( argv[3] );	
  if ( argv3 >> pavg){}
  else
    {
      cout<<"pavg should be a double"<<endl;
      return 1;
    }
  
  istringstream argv4( argv[4] );
  if ( argv4 >> range){}
  else
    {
      cout<<"range should be a double"<<endl;
      return 1;
    }
 
	     
  istringstream argv5( argv[5] );
  if ( argv5 >> max_failed_cws){}
  else
    {
      cout<<"max_failed_cws should be an int"<<endl;
      return 1;
    }

  //cout<<max_failed_cws<<endl;
  if (max_failed_cws<0)
    {
      max_num_of_cws=-max_failed_cws;
      max_failed_cws=1000000;
    }
 

  istringstream argv6( argv[6] );
  if ( argv6 >> lmax){}
  else
    {
      cout<<"lmax should be an int"<<endl;
      return 1;
    }

  
  istringstream argv8( argv[8] );
  if ( argv8 >> debug){}
  else
    {
      cout<<"debug should be an int"<<endl;
      return 1;
    }

  istringstream argv9( argv[9] );
  if ( argv9 >> channel){}
  else
    {
      cout<<"channel should be an int"<<endl;
      return 1;
    }

   istringstream argv10( argv[10] );
  if ( argv10 >> alpha){}
  else
    {
      cout<<"alpha should be double"<<endl;
      return 1;
    }
 istringstream argv11( argv[11] );
   if ( argv11 >> decode_p){}
  else
    {
      cout<<"decode_p should be double"<<endl;
      return 1;
    }
 istringstream argv12( argv[12] );
    if ( argv12 >> decode_prange){}
  else
    {
      cout<<"decode_prange should be double"<<endl;
      return 1;
    }

     istringstream argv13( argv[13] );
    if ( argv13 >> lambda){}
  else
    {
      cout<<"lambda should be double"<<endl;
      return 1;
    }

  double num_iter=0.0; //for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
  // int num_of_x_suc_dec=0;//number of Hx successfully decoded results
  int max_fail=0;//number of fails that reach maximum iterations
  int syn_fail=0;//number of fails that get the right syndrome but fails
  
  bool Hx_suc=false;
  // bool Hz_suc=false;



  //read the parity check matrix:
  int n1,n2,n,k,r1,r2,r;
  
  
  GF2mat Hx=read_matrix( n1,r1, file_name);
  GF2mat Hz=read_matrix( n2,r2, file_name2);
  //Hx=Hx.get_submatrix(0,0,r1-1,n1-1);
  //Hz=Hz.get_submatrix(0,0,r2-1,n2-1);
  //are the parity check matrices right?
  if (n1!=n2){cout<<"nx!=nz, the two matrices donot match"<<endl;return 1;}  
  n=n1;
  
   if (pavg>=1)
    {
      wt=pavg;
      pavg=1.0*wt/n;
     
      cout<<"pmin=wt/n*"<<1-range <<"pmax=wt/n*"<<1+range<<":"<<endl;  
      cout<<"for wt "<<wt<<" errors:"<<endl;
    }
   else
     {
       wt=0;// for p_min, pmax decode
     }
   if (range<1)
     {
  pmin=1.0*pavg*(1-range);
     }
   else {pmin=0;}
   
  pmax=1.0*pavg*(1+range);
  if (decode_prange<1)
    {
  decode_pmin=1.0*decode_p*(1-decode_prange);
    }

  else{decode_prange=0;}
  decode_pmax=1.0*decode_p*(1+decode_prange);
  r=r1+r2;
  int rankx=GF2mat_rank(Hx);
  int rankz=GF2mat_rank(Hz);
  k=n-rankx-rankz;
  GF2mat zero_mat1(r1,r2);

 
  

  if((Hx*Hz.transpose())==zero_mat1) {}
  else{cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl;return 1;}
  //cout<<Hx*Hz.transpose()<<endl;
  
  GF2mat Gx=get_gen(Hx);
  //cout<<Gx*(Hx.transpose())<<endl;
  GF2mat Gz=get_gen(Hz);
  // cout<<Gz*(Hz.transpose())<<endl;
  
  nodes  xchecks[r1];//checks for Hx and Z errors
  nodes  zerrors[n];
  nodes  zchecks[r2];//checks for Hz and X errors
  nodes  xerrors[n];

  int E1=0;
  int E2=0;
  //find the neighbourhoods of all nodes:
  initialize_checks (Hx, xchecks,  E1);
  initialize_errors(Hx, zerrors);

  initialize_checks (Hz, zchecks,  E2);
  initialize_errors(Hz, xerrors);

  vec px(n);
  vec px_dec(n);
  vec LR(n);
  int OSD_suc=0;
  // int  num_of_suc_dec=0;
  
         
  pro_dist( pmin,pmax, px);
  pro_dist( decode_pmin,decode_pmax, px_dec);  
  
  // int er=0;  //er is the number of one-bit errors (x for 1, z for 1, y for 2) that are wrong after decoding 
  while (  num_of_failed_cws<max_failed_cws&&num_of_cws<max_num_of_cws)
    {
      //cout<<num_of_cws<<endl;
      //cout<<max_num_of_cws<<endl;
      num_of_cws++;
      Hx_suc= quan_decode(Hx,Gz, xchecks,zerrors,px,px_dec,pmin,pmax,num_iter,lmax,wt,max_fail,syn_fail,debug,LR,rankx,OSD_suc,alpha,lambda);
 
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

   

   
  if ((debug/2)%2==1)
    {
      cout<<"parallel decoding"<<endl;
    }
  else
    {
      cout<<"serial decodeing"<<endl;
    }

  if ((debug/8)%2==1)
    {
      cout<<"use OSD if failed:"<<endl;
    }
  
  cout<<"lambda="<<lambda<<endl;
  cout<<"alpha="<<alpha<<endl;
  cout<<"real_ p=  ( "<<pmin<<", "<<pmax<<"),decode_p=("<<decode_pmin<<", "<<decode_pmax<<"), there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a [["<<n<<", "<<k<<"]] code (decode x errors only)"<<endl;
   cout<<"average iterations:"<<endl;
 
   cout<<num_iter/num_of_suc_dec<<"\n\n"<<endl;
   cout<<"syn_fail="<<syn_fail<<endl;
   cout<<"max_fail="<<max_fail<<endl;
   cout<<"OSD_suc="<<OSD_suc<<endl;

  
   // cout<<"num of zero errors is about "<<pow(p,n)*num_of_cws<<endl;
 
 
      
   ofstream myfile;
   myfile.open (data_file,ios::app);
   if (wt==0)
     {
       myfile << n<<" "<<d<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<pavg<<" "<<range<<" "<<1.0*num_iter/num_of_suc_dec<<"  "<<num_of_suc_dec<<" "<<num_of_cws<<"  "<<syn_fail<<" "<<max_fail<<" "<<1.0*syn_fail/num_of_cws<<" "<<1.0*max_fail/num_of_cws<<" "<<decode_p<<"  "<<decode_prange<<" "<<alpha<<" "<<OSD_suc<<" "<<lambda<<endl;
     }
   else
     {
       myfile << n<<" "<<d<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<wt<<" "<<range<<" "<<1.0*num_iter/num_of_suc_dec<<"  "<<num_of_suc_dec<<" "<<num_of_cws<<"  "<<syn_fail<<" "<<max_fail<<" "<<1.0*syn_fail/num_of_cws<<" "<<1.0*max_fail/num_of_cws<<decode_p<<"  "<<decode_prange<<" "<<alpha<<" "<<OSD_suc<<" "<<lambda<<endl;
     }
   myfile.close();
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "\n Run-time " << time_span.count() << " seconds.\n";
     
  return 0;

}
  

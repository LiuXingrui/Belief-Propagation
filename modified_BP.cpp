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


#ifndef modified_BP
#define modified_BP
#include"modified_BP.h"
#endif


using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

//Perform the matrix reduction for H and D, for stabilizer codes, if H=Hx, D=Hz, I call it D because it represents the degeneracy
void Mat_Trans(GF2mat& H, GF2mat& D, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug){

  GF2mat input_H=H;
  GF2mat input_D=D;
  GF2mat output_H, output_D;
  for (int i=0;i<num_of_row_reduction;i++)
    {
      row_reduction(i,input_H, input_D, output_H,output_D,A);
      //  check if the prog is right:
      if (debug==1)
	{
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  if((output_H*output_D.transpose())==zero_mat1) {}
	  else{cout<<"(output_H*output_D^T!=0 for the row "<<i<<endl;return;}
	}
      
      input_H=output_H;
      input_D=output_D;
    }
  H_tilde=output_H;
  D_tilde=output_D;
}


void row_reduction(int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<vector<int>>& A, int debug){
  vector<int> Ai;
  vector<int> Bi;
  int n=D.cols();
  int k=D.rows();
  int w=0;

  for (int j=0;j<n;j++)
    {
      if (D(which_row,j)==1){Ai.push_back(j);w++;}  // find where are the 1s
    }

  Ai.push_back(w);
  A.push_back(Ai);
  if (w>1)
    {
      Add_cols(H,D,Ai,w,debug);// now H -> H*A^T^(-1), D-> DA
      if (w>2)
	{
	  Add_cols_and_rows(H,D,Ai,w,H2,D2,debug);  // now H-> \frac{H*A^T^(-1)}{F}, D-> DA | DAB
	}
}

void Add_cols (GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug){

  int k=D.rows();
  int r=H.rows();

  for (int i=1;i<w;i++)
    {
      for (int j=0;j<k;j++)
	{
	  D.set(j,Ai[i],D(j,Ai[i])+D(j,Ai[0]));
	}
    }
  for (int j=0;j<r;j++)
    {
      H.set(j,0,0);
    }
  
  if (debug==1)
    {
      GF2mat zero_mat1(H.rows(),D.rows());
      if((H*D.transpose())==zero_mat1) {}
      else{cout<<"(after add_cols, H*D^T!=0,"<<endl;return;}
    }
      
}


void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug){

  int num_of_added_cols=pow(2,w)-w-1;
  int r=D.cols();
  GF2mat B,F;
  construct_BF(Ai,w,c,B,F);
  GF2mat D2_right=D*B;
  D2=merge_mat_hori(D,B);
  H2=merge_mat_vert(H,F);

  if (debug==1)
    {
      GF2mat zero_mat1(H2.rows(),D2.rows());
      if((H2*D2.transpose())==zero_mat1) {}
      else{cout<<"(after add_cols_and_rows, H*D^T!=0,"<<endl;return;}
    }
      
}

GF2mat construct_BF(vector<int> Ai, int w,int c,GF2mat&B,GF2mat& F){
  GF2mat B(c,pow(2,w-1)-w);

  for (i=2;i<=w;i++)
    {
      for (j

  


}





GF2mat merge_mat_vert(const GF2mat &up,const GF2mat &bottom, int debug)
{
  if (up.cols()!=bottom.cols())
    {
      cout<<"the 2 matrices cannot merge vertically"<<endl;
      GF2mat error(1,1); 
      return error;
    }

  else
    {
      int c=up.cols();
      int r=up.rows()+bottom.rows();
      int r1=up.rows();
      int r2=bottom.rows();
      GF2mat m(r,c);

      for (int i=0;i<c;i++)
	{
	  for (int j1=0;j1<r1;j1++){
	    m.set(j1,i,up(j1,i));
	  }
	  for (int j2=0;j2<r2;j2++){
	    m.set(j2+r1,i,bottom(j2,i);
	  }	  
	}
	  if (debug==1){cout<<"up:\n"<<up<<"\n  bottom:\n"<<bottom<<"\n result: \n"<<m<<endl;
       return m;
    }
}


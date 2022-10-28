#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <bits/stdc++.h>

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
void Mat_Trans(GF2mat& H, GF2mat& D, vector<double> &K, vector<double>& K_tilde, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug){

  GF2mat input_H=H;
  GF2mat input_D=D;
  GF2mat output_H, output_D;
  vector<double> input_K=K;
  vector<double> output_K=K;
  for (int i=0;i<num_of_row_reduction;i++)
    {
      //perform row reduction for 0th column of D
      row_reduction(0,input_H, input_D, output_H,output_D,input_K,output_K,A);//I havn't delete the rows and columns
 
      //  check if the prog is right:
      if (debug==1)
	{
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  cout<<"after "<<i<<"th row reduction before deleting rows/cols H:\n"<<output_H<<" \n D: \n"<<output_D<<endl;
	  if((output_H*output_D.transpose())==zero_mat1) {}
	  else{cout<<"(output_H*output_D^T!=0 for the row "<<i<<endl;return;}
	}
      //delete 0th col and row of D, delete 0th col of H
      output_H=output_H.get_submatrix(0,1,output_H.rows()-1,output_H.cols()-1);
      output_D=output_D.get_submatrix(1,1,output_D.rows()-1,output_D.cols()-1);
      
      if (debug==1)
	{
	
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  if((output_H*output_D.transpose())==zero_mat1) {}
	  else{cout<<"(after deleting rows and cols: output_H*output_D^T!=0 for the row "<<i<<endl;return;}
	}

      
      input_H=output_H;
      input_D=output_D;
      input_K=output_K;
    }
  H_tilde=output_H;
  D_tilde=output_D;
  K_tilde=output_K;
}


void row_reduction(int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<double> &input_K,vector<double>&output_K,vector<vector<int>>& A, int debug){
  vector<int> Ai;
  vector<vector<int>> b;
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
	  Add_cols_and_rows(H,D,Ai,w,H2,D2,debug,b);  // now H-> \frac{H*A^T^(-1)}{F}, D-> DA | DAB
	}
      K_trans(input_K,output_K,Ai,w,b);
    }
}
  // b stores the binary form of 1,2...2^(w-1)-1
void K_trans(vector<double>& input_K, vector<double> &output_K, vector<int>& Ai, int w,vector<vector<int>>& b){
  int num_change_cols=pow(2,w-1)-1;
  int num_B_tau=pow(2,w);
  int change_col=0;
  vector<double> B_tau;
 
  vector<vector<int>>  Tau;
 
  // get the binray form of tau
   for (int i=0;i<num_B_tau;i++)
   { 
      vector<int> tau;  
      for (int j=0;j<w;j++)
	{
	  int twotoj=pow(2,j)+0.5;
	  //temp is jth digit of i
	 int temp=i%twotoj;
	  // if jth digit is one, set corresponding entries in co and ro to be one
	  if (temp==1){tau.push_back(1);}
	  else {tau.push_back(0);};
	}
	Tau.push_back(tau);   
   }   
  construct_B_tau(B_tau,w,input_K,Ai,Tau);   
	
  for (int i=1;i<=num_change_cols;i++) // every changed col gives a changed K
    {
      double temp=1;
      vector<int> bi=b[i-1]; // bi is the subscript of the current K
      int wt_bi=bi[bi.size()-1];
      vector<int> bi_to_tau_pos;
	  
      //because the vector bi has weight w-1, but the vector tau has weight w, every even weight tau gives a bi, 
     // when convert bi back to tau,just a a digit at the end of bi, if bi has even weight, the first digit of tau=0, if bi has odd weight, the first digit of tau=1.
      if (wt_bi%2==0){}
      else {bi_to_tau_pos.push_back(0);}  // if bi has odd weight, the first digit of tau is 1.
      
      for (int k=0;k<bi.size()-1;k++)
	{
	  if (bi[k]==1){bi_to_tau_pos.push_back(k+1);}  //find the entries of tau with 1
	}
      
      for (int j=0;j<num_B_tau;j++)
	{
	  int wt_tau_bi=0;
	  for (int ii=0;ii<bi_to_tau_pos.size();ii++)
	    {
	      if (Tau[j][bi_to_tau_pos[ii]]==1){wt_tau_bi++;}//Tau[j] is the binary form of j, check the number of i such that Tau[j][i]=tau[i]=1
	    }
	      if (wt_tau_bi%2==0){temp=temp*B_tau[j];}
	      else{temp=temp/B_tau[j];}
	}      
      if (bi[bi.size()-1]==1){output_K[Ai[change_col]]=1/pow(2,w)*log(temp);change_col++;}  //change the value of corresponding K
      else {output_K.push_back(1/pow(2,w)*log(temp));}
    }  
}
      
void construct_B_tau(vector<double> &B_tau,int w, vector<double>& input_K, vector<int> Ai,vector<vector<int>>&  Tau){
   int num_B_tau=pow(2,w);
   int temp=0;
   for (int i=0;i<num_B_tau;i++)
     {
       temp=0;
       for (int j=0;j<w;j++)
	 {
	   temp=temp+input_K[Ai[j]]*pow(-1,Tau[i][j]);
	 }     
       B_tau.push_back(temp);
     }	 	 
 }
 
void Add_cols(GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug){

  int k=D.rows();
  int r=H.rows();

  for (int i=1;i<w;i++)  // choose 1st, 2nd to w-1 cols, and 0th column to them,
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


void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug,vector<vector<int>>& b){

  int num_add_cols=pow(2,w-1)-w; //number of added cols, also is the num of cols of B,cols of B are binray form of 0,1,...2^(w-1)-1,
                                  //but delete w-1 weight=1 cols, and 1 weight=0 col, so num_add_col=2^(w-1)-(w-1)-1=2^(w-1)-w
  int c=D.cols();
  int r=D.rows();
  GF2mat  DB(c,num_add_cols);
  GF2mat  F(num_add_cols,H.cols()+num_add_cols);
  
  construct_BF(D,Ai,w,c,DB,F,b);
  
  GF2mat zero_mat(H.rows(),num_add_cols);
  GF2mat H1=merge_mat_hori(H,zero_mat);
  
  D2=merge_mat_hori(D,DB);
  H2=merge_mat_vert(H1,F);

  if (debug==1)
    {
      GF2mat zero_mat1(H2.rows(),D2.rows());
      if((H2*D2.transpose())==zero_mat1) {}
      else{cout<<"(after add_cols_and_rows, H*D^T!=0,"<<endl;return;}
    }
      
}

// b stores the binary form of 0...2^(w-1)-1, and the last element is w
void construct_BF(GF2mat& D,vector<int> &Ai, int w,int c, GF2mat& DB, GF2mat& F,vector<vector<int>>& b){

  int added_c=0;
  int r=D.rows();
  int num_add_cols=pow(2,w-1)-w;
  int num_change_cols=pow(2,w-1)-1; // number of changed cols,which equals to num_add_cols_+w-1
  GF2mat B(c,num_add_cols);
  bvec co(c); // the column added to the right of D
  bvec ro(F.cols());//the row added to the bottom of H
  ro.zeros();
  co.zeros();
  int col_ind=0;  // a column index
  int temp=0;

  
  
  //for B:
  // i starts with 1, because i=0 gives the all-zero col.
  for (int i=1;i<=num_change_cols;i++)
    {
      int wt=0;
      co.zeros();
      ro.zeros();
      vector<int> bi;  //bi is the binary form of i
      //calculate bi, and where are the corresponding 1-entries of co
      for (int j=1;j<=w-1;j++)
	{
	  //temp is jth digit of i
	  int twotoj=pow(2,j-1)+0.5;
	  temp=i%twotoj;
	  // if jth digit is one, set corresponding entries in co and ro to be one
	  if (temp==1){ co(Ai[j-1])=1;ro(Ai[j-1])=1;bi.push_back(1);wt++;}
	  else {bi.push_back(0);};
	}

      //if wt=1, we don't need to add this col
      if (wt>1)
	{
	  B.set_col(col_ind,co);
	  col_ind++;
	}
      
      if (wt==2){F.set_row(col_ind,ro);}

      // if wt>2, we can set F(i,i)=F(i,i-1)=1 first and to find  which entry should equals to 1 that keeps the orthogonality.
      if (wt>2)
	{
	  F.set(col_ind,col_ind-1,1);
	  bvec temp_vec1=B.get_col(col_ind);
	  bvec temp_vec2=B.get_col(col_ind-1);
	  
	  // substract the corrsponding columns of B to find where is the third col.
	  bvec temp_vec21=temp_vec1+temp_vec2;
	  for (int k=0;k<col_ind-1;k++)
	    {
	      if (temp_vec21==B.get_col(k)){F.set(col_ind,k,1);break;}
	      if (k==col_ind-2){cout<<"k==col_ind-2, something went wrong"<<endl;return;}
	    }
	}
          
      F.set(col_ind,col_ind,1);
      bi.push_back(w);// store the weight in the last position of bi
      b.push_back(bi);
    }
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
	  for (int j1=0;j1<r1;j1++)
	    {
	      m.set(j1,i,up(j1,i));
	    }
	  for (int j2=0;j2<r2;j2++)
	    {
	      m.set(j2+r1,i,bottom(j2,i));
	    }	  
	}
      if (debug==1){cout<<"up:\n"<<up<<"\n  bottom:\n"<<bottom<<"\n result: \n"<<m<<endl;}
      return m;
    }
}

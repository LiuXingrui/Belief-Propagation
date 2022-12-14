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
void Mat_Trans(const GF2mat& H, const GF2mat& D, vector<double> &K, vector<double>& K_tilde, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug){

  GF2mat input_H=H;
  GF2mat input_D=D;
  GF2mat output_H, output_D;
  vector<double> input_K=K;
  vector<double> output_K=K;
  // cout<<"start"<<endl;
  
  int real_num_of_row_reduction;
  if (num_of_row_reduction>=0){  real_num_of_row_reduction=num_of_row_reduction;}
  else{real_num_of_row_reduction=-num_of_row_reduction;} //negative num_of_row_reduction for combine e at the bottom of D


    if (debug==1)
	{
	
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  cout<<"before row reduction H is"<<H<<" \n D: \n"<<D<<endl;
	  
	  cout<<"K_tilde is \n:"<<endl;
	  for(auto ii:input_K)
	    {
	      cout<<ii<<" ";
	    }
	  cout<<endl;
	}
      
  for (int i=0;i<real_num_of_row_reduction;i++)
    {
      // cout<<i<<"th row reduction"<<endl;
      int min_wt=H.cols();
      int this_row=0;

      if (num_of_row_reduction>=0)
	{
	  for (int j=0;j<input_D.rows();j++)
	    {
	  
	      bvec temp_vec=input_D.get_row(j);
	      int temp_wt=Weight(temp_vec);
	      if (temp_wt<min_wt){min_wt=temp_wt;this_row=j;}
	    }
	}
      
      else
	{
	  for (int j=0;j<input_D.rows()-1;j++) //donot touch the last row- the error vector
	    {
	  
	      bvec temp_vec=input_D.get_row(j);
	      int temp_wt=Weight(temp_vec);
	      if (temp_wt<min_wt){min_wt=temp_wt;this_row=j;}
	    }
	}
      
   
      if (debug==1)
	{
	  cout<<i+1<<"th row reduction, delete "<<this_row<<"th row:"<<endl;
	}
      //perform row reduction for the column with minimum weight of D
      // cout<<121<<endl;
      if (min_wt>15){cout<<"for "<<i<<"th row reduction (first index is 0) weight of minimum-weight row of D is "<<min_wt<<", so stop here"<<endl;break;}
      row_reduction(this_row,input_H, input_D, output_H,output_D,input_K,output_K,A,debug);//I havn't delete the rows and columns
      
      //  cout<<22<<endl;
      //  check if the prog is right:
      if (debug==1)
	{
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  //cout<<"after "<<i<<"th row reduction before deleting rows/cols \n H:\n"<<output_H<<" \n D: \n"<<output_D<<endl;
	  if((output_H*output_D.transpose())==zero_mat1) {}
	  else{cout<<"before deleting: output_H*output_D^T!=0 for the row "<<i<<endl;}
	}

      

      if (A[i][0]==0) //if delete the 0th column from H and D
	{
	  output_H=output_H.get_submatrix(0,1,output_H.rows()-1,output_H.cols()-1);
	  if (i!=D.rows()-1) { output_D=output_D.get_submatrix(0,1,output_D.rows()-1,output_D.cols()-1);      }
	}
      
      else if(A[i][0]==output_H.cols()-1)//if delete the last column
	{
	  output_H=output_H.get_submatrix(0,0,output_H.rows()-1,output_H.cols()-2);
	  if (i!=D.rows()-1) { output_D=output_D.get_submatrix(0,0,output_D.rows()-1,output_D.cols()-2);}
	}

      
      else //delete A[i][0]th column from H and D
	{
	  GF2mat H1=output_H.get_submatrix(0,0,output_H.rows()-1,A[i][0]-1);
	  GF2mat H2=output_H.get_submatrix(0,A[i][0]+1,output_H.rows()-1,output_H.cols()-1);
	  output_H=merge_mat_hori(H1,H2);
	  if (i!=D.rows()-1)
	    {
	      GF2mat D1=output_D.get_submatrix(0,0,output_D.rows()-1,A[i][0]-1);
	      GF2mat D2=output_D.get_submatrix(0,A[i][0]+1,output_D.rows()-1,output_D.cols()-1);
	      output_D=merge_mat_hori(D1,D2);
	    }	  
	}
      if (i!=D.rows()-1)
	{
	  if (this_row==0) //if delete 0th row from D
	    {
	      output_D=output_D.get_submatrix(1,0,output_D.rows()-1,output_D.cols()-1);
	    }
	  else if(this_row==output_D.rows()-1)
	    {
	      output_D=output_D.get_submatrix(0,0,output_D.rows()-2,output_D.cols()-1);
	    }
	  else 
	    {
	      GF2mat D3=output_D.get_submatrix(0,0,this_row-1,output_D.cols()-1);
	      // cout<<D3.cols()<<endl;
	      GF2mat D4=output_D.get_submatrix(this_row+1,0,output_D.rows()-1,output_D.cols()-1);
	      // cout<<D4.cols()<<endl;
	      output_D=merge_mat_vert(D3,D4);
	    }
     
	}
      output_K.erase(output_K.begin()+this_row);

      //need delete an element from K
      
      if (debug==1)
	{
	
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  cout<<"after deleting "<<i+1<<" rows (delete row: "<<this_row<<" and column: "<<A[i][0]<<" )\n H:\n"<<output_H<<" \n D: \n"<<output_D<<endl;
	  if (i!=D.rows()-1)
	    {
	      if((output_H*output_D.transpose())==zero_mat1) {}
	      else{cout<<"after deleting rows and cols: output_H*output_D^T!=0 for the row the result is"<<output_H*output_D.transpose()<<endl;return;}
	    }
	  
	  cout<<"K_tilde is \n:"<<endl;
	  for(auto ii:output_K)
	    {
	      cout<<ii<<" ";
	    }
	  cout<<endl;
	}
      
    
      
      input_H=output_H;
      input_D=output_D;
      input_K=output_K;
    }

  H_tilde=input_H;
  D_tilde=input_D;
  K_tilde=input_K;
  // cout<<"end"<<endl;
}


void row_reduction(const int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<double> &input_K,vector<double>&output_K,vector<vector<int>>& A, int debug){
  vector<int> Ai;
  vector<vector<int>> b;
  int n=D.cols();
  int k=D.rows();
  int w=0;

  for (int j=0;j<n;j++)
    {
      if (D(which_row,j)==1){Ai.push_back(j);w++;}  // find where are the 1s
    }

  //cout<<1<<"wt is "<<w<<endl;
  // if (w==72){cout<<D<<endl;}
  Ai.push_back(w);
  A.push_back(Ai);
  if (w>1)
    {
      Add_cols(H,D,Ai,w,debug);// now H -> H*A^T^(-1), D-> DA
      //  cout<<2<<endl;
      if (w>2)
	{
	  // cout<<"before add cols and rows"<<endl;
	  Add_cols_and_rows(H,D,Ai,w,H2,D2,debug,b);  // now H-> \frac{H*A^T^(-1)}{F}, D-> DA | DAB
	  //  cout<<3<<endl;
	  // cout<<"after add cols and rows"<<endl;
	}
      else{H2=H;D2=D;}
      K_trans(input_K,output_K,Ai,w,b);
      //  cout<<4<<endl;
      // cout<<"after K_trans"<<endl;
    }
  else{H2=H;D2=D;}
}

  // b stores the binary form of 1,2...2^(w-1)-1
void K_trans(vector<double>& input_K, vector<double> &output_K, vector<int>& Ai, int w,vector<vector<int>>& b){
  int num_change_cols=pow(2,w-1)-1;// w=2: 1; w=3: 3; w=4:7
  int num_B_tau=pow(2,w);
  int change_col=0;
  vector<double> B_tau;
 
  vector<vector<int>>  Tau;

  int temp_K_changed=0;
  
  // get the binray form of i, store in the vector called "tau",notice the ith element in "tau" is the ith digit of i, so if i=1010, tau=[0,1,0,1]
  for (int i=0;i<num_B_tau;i++)
    {
      int tempi=i;
      vector<int> tau;  
      for (int j=0;j<w;j++)
	{
	  int tempi_red=tempi%2;
	  tempi=tempi/2;
	  
	  // if jth digit is one, set corresponding entries in co and ro to be one
	  if (tempi_red==1){tau.push_back(1);}
	  else {tau.push_back(0);};
	}
      Tau.push_back(tau);   
    }
  //now construct all B_tau
  construct_B_tau(B_tau,w,input_K,Ai,Tau);
  /*
  cout<<"w is "<<w<<"size of B_tau is "<<B_tau.size()<<endl;
  cout<<"B_tau is"<<endl;
  
  for (auto ii:B_tau)
    {
      cout<<ii<<" ";
    }
  cout<<endl;
  */
  
  for (int i=1;i<=num_change_cols;i++) // every changed col gives a changed K
    {
      double tempK=1;
      vector<int> i_binary_pos_plus_1; //so that is the position of corresponding tau for K
      int wt_i=0;
      int tempi=i;
      //get the binary form of i,store the position of 1 in "i_binary_pos"
      for (int ii=0;ii<=w-2;ii++)
	{
	  int tempi_red=tempi%2;
	  tempi=tempi/2;
	  if (tempi_red==1){i_binary_pos_plus_1.push_back(ii+1);wt_i++;}
	  
	}
      //get the tau for K
      vector<int> tau_for_K_pos=i_binary_pos_plus_1;
      if (wt_i%2==1) {tau_for_K_pos.push_back(0);}

	  
	  
      //because the vector i has weight w-1, but the vector tau has weight w, every even weight tau gives a i, 
      // when convert i back to tau,just add a digit at the end of i, if i has even weight, the first digit of tau=0, if bi has odd weight, the first digit of tau=1.
      double temp=1;
      for (int j=0;j<num_B_tau;j++)
	{
	  int wt_tau_bi=0;
	  for (int ii=0;ii<tau_for_K_pos.size();ii++)
	    {
	      if (Tau[j][tau_for_K_pos[ii]]==1){wt_tau_bi++;}//Tau[j] is the binary form of j, check the number of i such that Tau[j][i]=tau[i]=1
	    }
	  if (wt_tau_bi%2==0){temp=temp*B_tau[j];}
	  else{temp=temp/B_tau[j];}
	}      
      if (wt_i==1){output_K[Ai[temp_K_changed]]=1/pow(2,w)*log(temp);temp_K_changed++;}  //change the value of corresponding K
	else {output_K.push_back(1/pow(2,w)*log(temp));}
      }
    
  
}
      
void construct_B_tau(vector<double> &B_tau,int w, vector<double>& input_K, vector<int> Ai,vector<vector<int>>&  Tau){
   int num_B_tau=pow(2,w);
   double temp=0;
   for (int i=0;i<num_B_tau;i++)
     {
       temp=0;
       for (int j=0;j<w;j++)
	 {
	   temp=temp+input_K[Ai[j]]*pow(-1,Tau[i][j]);
	 }
       temp=cosh(temp);
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
      H.set(j,Ai[0],0);
    }
  
  if (debug==1)
    {
      GF2mat zero_mat1(H.rows(),D.rows());
      cout<<"after add cols \n H:\n"<<H<<" \n D: \n"<<D<<endl;
      if((H*D.transpose())==zero_mat1) {}
      else{cout<<"after add_cols, H*D^T!=0,"<<endl;return;}
    }
      
}


void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug,vector<vector<int>>& b){

  int num_add_cols=pow(2,w-1)-w; //number of added cols, also is the num of cols of B,cols of B are binray form of 0,1,...2^(w-1)-1,
                             //but delete w-1 weight=1 cols, and 1 weight=0 col, so num_add_col=2^(w-1)-(w-1)-1=2^(w-1)-w
  int c=D.cols();
  int r=D.rows();
  int n=H.cols();
  GF2mat  DB(c,num_add_cols);
  GF2mat  F(num_add_cols,H.cols()+num_add_cols);
  
  construct_BF(D,Ai,w,c,DB,F,b,n,debug);
  
  GF2mat zero_mat(H.rows(),num_add_cols);
  GF2mat H1=merge_mat_hori(H,zero_mat);


  D2=merge_mat_hori(D,DB);
  H2=merge_mat_vert(H1,F);
  
  if (debug==1)
    {
      // GF2mat zero_mat1(H2.rows(),D2.rows());
      //if((H2*D2.transpose())==zero_mat1) {}
      // else{cout<<"after add_cols_and_rows, H*D^T!=0, that is \n   H\n"<<H2<<"\n D: \n"<<D2<<"\n  result:\n"<<H2*D2.transpose()<<endl;return;}
    }
      
}

// b stores the binary form of 0...2^(w-1)-1, and the last element is w
void construct_BF(GF2mat& D,vector<int> &Ai, int w,int c, GF2mat& DB, GF2mat& F,vector<vector<int>>& b,int n,int debug){

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
      int ii=i;
      co.zeros();
      ro.zeros();
      vector<int> bi;  //bi is the binary form of i
      //calculate bi, and where are the corresponding 1-entries of co
      for (int j=1;j<=w-1;j++)
	{
	  //temp is jth digit of i
	  temp=ii%2;
	  ii=ii/2;
	  // cout<<"j is "<<j<<endl;
	  //   cout<<"temp is "<<temp<<endl;
	  //  cout<<"ii is "<<ii<<endl;
	
	  // if jth digit is one, set corresponding entries in co and ro to be one
	  if (temp==1){ co(Ai[j])=1;ro(Ai[j])=1;bi.push_back(1);wt++;}
	  else {bi.push_back(0);};
	}

      //if wt=1, we don't need to add this col
      if (wt>1)
	{
	  //cout<<"col_ind is"<<col_ind<<endl;
	  B.set_col(col_ind,co);
	  if (wt==2){F.set_row(col_ind,ro);}

	  // if wt>2, we can set F(i,i)=F(i,i-1)=1 first and to find  which entry should equals to 1 that keeps the orthogonality.
	  if (wt>2)
	    {
	     
	      //cout<<"B is\n"<<B<<endl;
	      // cout<<"F is\n"<<F<<endl;
	      F.set(col_ind,n+col_ind-1,1);
	      
	      bvec temp_vec1=B.get_col(col_ind);
	      bvec temp_vec2=B.get_col(col_ind-1);
	  
	      // substract the corrsponding columns of B to find where is the third col.
	      bvec temp_vec21=temp_vec1+temp_vec2;
	      int temp_wt=0;
	      int temp_one;
	      //if wt of temp_vec21 =1, then temp_vec1 is the sum of temp_vec2 and a column from D
	      //if wt >1, is the sum of temp_vec2 and a column from DB
	      for (int iii=1;iii<w;iii++)
		{
		  if (temp_vec21(Ai[iii])==1){temp_wt++;temp_one=iii;}
		}
	      if (temp_wt==1){F.set(col_ind,Ai[temp_one],1);}
	      else
		{
		  for (int k=0;k<col_ind-1;k++)
		    {
		      if (temp_vec21==B.get_col(k)){F.set(col_ind,k+n,1);break;}
		      if (k==col_ind-2){cout<<"k==col_ind-2, something went wrong"<<endl;return;}
		    }
		}
	    }
	  //cout<<3333<<endl;
	  F.set(col_ind,n+col_ind,1);
	  col_ind++;
	}
      bi.push_back(w);// store the weight in the last position of bi
      b.push_back(bi);
    }
  DB=D*B;

  
  //if (debug==1){
  // cout<<"B is \n"<<B<<endl;
  // cout<<"F is\n"<<F<<endl;
  // }
  
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
      //if (debug==1){cout<<"up:\n"<<up<<"\n  bottom:\n"<<bottom<<"\n result: \n"<<m<<endl;}
      return m;
    }
}


bool new_decoder1(GF2mat& H2,GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,double& num_iter, int lmax,int& max_fail,int&syn_fail,int debug, vec &LR,int rankH,int options)
{
  int v=H.cols();
  int c=H.rows();

  int wt_real_e=0;
  GF2mat real_e(v,1);
  
  wt_real_e=error_channel(real_e, pv);
  if (wt_real_e==0)
    {	  
      return true;
    }
    
  //if no error, break

  GF2mat zero_rvec(1,v);
   
  LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat zero_rvec2(v-rankH,1);
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix

      GF2mat output_e(v,1);
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
   	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,1);	
	    
	   if (H*output_e==syndrome)
		{		  
		  if(G*(output_e+real_e)==zero_rvec2)
		    {
		      num_iter= num_iter+l;		  
		      return true;
		    }
		  else
		    {
		      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
		      syn_fail++;		
		      return false;      	        
		    }	    	  
		}  
	}
      /*
      cout<<"output e is \n"<<endl;
      err_pos2(output_e);
      cout<<"real e is \n"<<endl;
      err_pos2(real_e);      
      */      
      max_fail++;
      return false;
 }


bool new_decoder2(GF2mat& H,GF2mat& G,GF2mat &H_tilde, GF2mat &H_tilde_star,const vector<vector<int>>& A,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,const int&  num_row_red,double& num_iter, int lmax,int& max_fail,int&syn_fail,int& trans_suc,int& no_error,int debug, vec &LR,int rankH,int options)
{
  int v=H.cols();
  int c=H.rows();
  int vt=H_tilde.cols();
  int ct=H_tilde.rows();
  if (num_row_red==0){if(H==H_tilde){}else {cout<<"wrong H_tilde"<<endl;}}

  int wt_real_e=0;
  GF2mat real_e(v,1);
  
  wt_real_e=error_channel(real_e, pv);
  if (wt_real_e==0)
    {
      trans_suc++;
      no_error++;
      return true;
    }
    
  //if no error, break

  GF2mat zero_rvec(1,v);
   
  LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat zero_rvec2(v-rankH,1);
  // cout<<11<<endl;
  GF2mat zero_vec_for_syndrome(ct-c,1);
  // cout<<12<<endl;
  
  GF2mat syndrome=H*real_e;


  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      trans_suc++;
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
  GF2mat syndrome_tilde(ct,1);
  // cout<<12.125<<endl;
  // cout<<ct-c<<endl;
  if (ct-c!=0){ syndrome_tilde=merge_mat_vert(syndrome,zero_vec_for_syndrome);}
  else{syndrome_tilde=syndrome;}
  // cout<<12.25<<endl;
  
  mat mcv(ct,vt);   // the messages from check nodes to variable nodes, which are p0/p1
  mat mvc(ct,vt);   // the messages from variable nodes to check nodes
  mcv.zeros();
  mvc.zeros();
  
  //  cout<<12.5<<endl;
  initialize_massages( mcv,mvc, H_tilde); //initialize to all-1 matrix
  
  //  cout<<13<<endl;
  
  GF2mat output_e(v,1);
  GF2mat output_e_tilde(vt,1);
      
  for (int l=1;l<=lmax;l++)
    {
      //  cout<<l<<endl;
      quan_s_update(checks,errors, mcv,mvc,syndrome_tilde,pv_dec, ct, vt,output_e_tilde,LR,1);
      
      //  cout<<14<<endl;
      
      if (H_tilde*output_e_tilde==syndrome_tilde)
	{
	  inverse_trans(output_e,output_e_tilde,A);
	  /*
	  if(num_row_red==0){if (output_e==output_e_tilde){}else{cout<<"wrong output_e"<<endl;}}
	  if(num_row_red==0){if (syndrome==syndrome_tilde){}else{cout<<"wrong syndrome"<<endl;}}
	  if(num_row_red==0){cout<<"H*output_e is \n"<<H*output_e<<"\n H*real_e is \n"<<H*real_e<<endl;}
	  */
	  trans_suc++;
	  //cout<<"ok"<<endl;
	  // return true;
	  if(H*(output_e+real_e)==zero_mat1)
	    {
	      num_iter= num_iter+l;		  
	      return true;
	    }
	  else
	    {
	      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
	      syn_fail++;		
	      return false;      	        
	    }	    	  
	}  
    }
	            
      max_fail++;
      return false;
 }


int Weight(const bvec &cw)
{
  int n=cw.size();
  int wt=0;
  for (int i=0;i<n;i++)
    {
      if(cw(i)==1){wt++;}
    }
  return wt;
}

void inverse_trans(GF2mat& output_e,const GF2mat& output_e_tilde,const vector<vector<int>>& A)
{
  // cout<<1<<endl;
  int num_row_red=A.size();  //
  int v=output_e.rows();
  int temp=0;
  vector<int> temp_vec;
  for (int j=0;j<output_e_tilde.rows();j++){temp_vec.push_back(output_e_tilde.get(j,0));}
  //  cout<<2<<endl;

  //here I give an example, assume G has a row (1111), we have an error e=(1011) or (0100) (they are degenerate, since H must have even number of 1s in this row),
  //after a row reduction, e-> (100) no matter which one it is. To convert it back, just add a zero at the position a column was deleted.
  
  for (int i=num_row_red-1;i>=0;i--)
    {
      temp=A[i][0];
      temp_vec.insert(temp_vec.begin()+temp,0);      
    }
  // cout<<3<<endl;
  for (int i=0;i<v;i++){output_e.set(i,0,temp_vec[i]);}

  // cout<<"inverse_trans: num_row_red="<<num_row_red<<endl;
  // if (num_row_red==0){cout<<"output_e is \n"<<output_e<<"\n e_tilde is \n"<<output_e_tilde<<endl;}
  
  //  cout<<4<<endl;
}

//here eT is a row vector, because in other functions usually error is a column vector
// something wrong here, so empty_G should always be true
double energy (bool empty_G,const GF2mat& G, const GF2mat& eT, const vector<double> &K)
{

  if (eT.rows()!=1){cout<<"eT.rows!=1, some thing wrong with the energy function"<<endl;return 0;}
  if (empty_G==true)
    {
      // cout<<"e1"<<endl;
      double temp=0;
      for (int i=0;i<eT.cols();i++)
	{
	  int temp2=eT.get(0,i);//becaues get return a bin, convert it to int
	  temp=K[i]*pow(-1,temp2)+temp;
	}
      //cout<<"e2"<<endl;
      return temp;
    }

  
  if(eT.cols()!=G.cols()){cout<<"error in energy func:sizes of e and G do not match"<<endl;return 0;}

  double temp=0;
  GF2mat tempmat=eT*G;
  for (int i=0;i<eT.cols();i++)
    {
      int temp2=tempmat.get(0,i);
      temp=K[i]*pow(-1,temp2)+temp;
    }
  return temp;
}

double ML_suc_rate (vec p,const GF2mat& D,const GF2mat& H, const GF2mat& H_tilde_star, vector<double>&K,int d,int num_large_wt_error)
{
  cout<<"ML_suc_rate start"<<endl;
  
  double suc_rate=0;
  int v=H.cols();
  int c=H.rows();
  
  int vt=H_tilde_star.cols();

  int wt_real_e=0;
  GF2mat real_e(v,1);
  int large_wt_error=0;
  int num_error=0;
  vector<double> K_tilde;

  while (large_wt_error< num_large_wt_error)
    {
      K_tilde.clear();
      num_error++;
      
      int wt_real_e=0;
      GF2mat real_e(v,1); 
      wt_real_e=error_channel(real_e, p);

      // cout<<0<<endl;
      
      if (wt_real_e<d/2){suc_rate++;}   
      else
	{
	  large_wt_error++;
	  GF2mat real_eT=real_e.transpose();
	  
	  GF2mat De=merge_mat_vert(D,real_eT);
	  vector<vector<int>> A;
	  GF2mat H_tilde,e_tilde,empty_tempmat;
	  
	  
	  //  cout<<1<<endl;
	  // here e_tilde is the result transformed matrix De, so its transpose is the error vector transformed
	  Mat_Trans(H,De,K,  K_tilde, H_tilde, e_tilde, A,-D.rows());

	  //	  cout<<2<<endl;
	  if (e_tilde.cols()==H_tilde_star.cols() && e_tilde.rows()==H_tilde_star.rows())
	    {
	      double exp_energy_ratio=exp(energy (true,empty_tempmat, e_tilde, K_tilde))/(exp(energy (true,empty_tempmat, e_tilde, K_tilde))+exp(energy (true,empty_tempmat, e_tilde+H_tilde_star, K_tilde)));
	      //   double exp_e_plus_c_energy=0;
	      
	      if (exp_energy_ratio>0.5){suc_rate++;}
	    }
	  else {cout<<"error: e_tilde.cols()="<<e_tilde.cols()<<" \n H_tilde_star.cols()="<<H_tilde_star.cols()<<"\n  e_tilde.rows()="<<e_tilde.rows()<<"\n H_tilde_star.rows()="<<H_tilde_star.rows()<<endl;}
	}
    }
  cout<<" number of codewords tested: "<<num_error<<endl;

  /*
  cout<<"K_tilde in ML calculating is \n"<<endl;
  for (int ii=0;ii<K_tilde.size();ii++)
    {
      cout<<K_tilde[ii]<<endl;
    }
  */
  
  suc_rate=suc_rate/num_error;
  return suc_rate;

  
}

void ML_decoder_verify(vec p,const GF2mat& H, const GF2mat D,const GF2mat& H_star, const GF2mat D_star,const GF2mat H_tilde, const GF2mat H_tilde_star,vector<double>&K,int max_num_cws,double& after_trans_suc_rate,double&suc_rate,double &after_trans_theoric_suc_rate,double& min_wt_suc_rate,int debug)
{
  // get logical operator:
  GF2mat L;
  // cout<<"ML_suc_rate start"<<endl;
  get_logical(H, D,H_star, L);

  int min_wt_suc_count=0;
  int after_trans_suc_count=0;
  int suc_count=0;
  int after_trans_theoric_suc_count=0;

  int num_cws_count=0;
  
   // bool after_trans_suc=true;
   // bool suc=true;
  
  int v=H.cols();
  int c=H.rows();
  
  int vt=H_tilde_star.cols();

  vector<double> K_tilde;

  //generate a random error
  while (num_cws_count<max_num_cws)
    {
      num_cws_count++;
      K_tilde.clear();
      
      int wt_real_e=0;
      GF2mat real_e(v,1); 
      wt_real_e=error_channel(real_e, p);

      if (debug==1){ cout<<"real_e is \n"<<real_e<<endl;}
 

      // cout<<0<<endl;
      
      GF2mat real_eT=real_e.transpose();
      
      GF2mat De=merge_mat_vert(D,real_eT);
      vector<vector<int>> A;
      GF2mat H_tilde,e_tilde,empty_tempmat;
	    
      // cout<<3<<endl;
      //here e_tilde is a row vector so it is in fact the transposed one
      Mat_Trans(H,De,K,  K_tilde, H_tilde, e_tilde, A,-D.rows());
           if (debug==1)
	     {
	       cout<<"after trans real_e is \n"<<e_tilde<<endl;
	       cout<<"e_tilde+logical_operattor is\n"<<e_tilde+H_tilde_star<<endl;
	     }

      //  cout<<4<<endl;

      //theoric ML suc rate:
      if (e_tilde.cols()==H_tilde_star.cols() && e_tilde.rows()==H_tilde_star.rows())
	{
	  double exp_energy_ratio=exp(energy (true,empty_tempmat, e_tilde, K_tilde))/(exp(energy (true,empty_tempmat, e_tilde, K_tilde))+exp(energy (true,empty_tempmat, e_tilde+H_tilde_star, K_tilde)));
	  //   double exp_e_plus_c_energy=0;
	      
	  if (exp_energy_ratio>0.5)
	    {
	      after_trans_theoric_suc_count++;
	           if (debug==1)
		     {
		       cout<<"in theory: suc"<<endl;
		     }
	    }
	}
	  else {cout<<"error: e_tilde.cols()="<<e_tilde.cols()<<" \n H_tilde_star.cols()="<<H_tilde_star.cols()<<"\n  e_tilde.rows()="<<e_tilde.rows()<<"\n H_tilde_star.rows()="<<H_tilde_star.rows()<<endl;}


      //ML decoder:
      GF2mat syndrome=H*real_e;
      GF2mat algebraic_e(H.cols(),1);
      GF2mat ML_output_e(H.cols(),1);
      //  cout<<5<<endl;
      
      algebraic_decoder(H,syndrome,algebraic_e);
     
      
      // cout<<6<<endl;
      //  cout<<"H is \n"<<H<<" \n H_star is \n"<<H_star
      if (debug==1)  {cout<<"ML decoder:\n algebraic_e is\n"<<algebraic_e<<endl; }
      // cout<<"ML decoder:\n algebraic_e is\n"<<algebraic_e<<endl;
 
      // cout<<"algebraic_e +L is\n"<<algebraic_e+L.transpose()<<endl;

      int alge_e_wt=s_weight(algebraic_e);
      int alge_e_L_wt=s_weight(algebraic_e+L.transpose());
      GF2mat min_wt_e;
      
      if (alge_e_wt<=alge_e_L_wt) {min_wt_e=algebraic_e;}
      
      else {min_wt_e=algebraic_e+L.transpose();}
      
      GF2mat tempzeromat(D_star.rows(),1);
      if (H*min_wt_e==syndrome)
	{
	  if (D_star*(min_wt_e+real_e)==tempzeromat)
	    {
	      min_wt_suc_count++;
	    }
	}
      
	  ML_decoder(algebraic_e, D,L,K,ML_output_e);
      
           if (debug==1)
	     {
	       cout<<"ML decoder: output_e is\n"<<ML_output_e<<endl;
	     }

      //  cout<<7<<endl;
    
      if (H*ML_output_e==syndrome)
	{
	  //  cout<<"D_star is\n"<<D_star<<endl;
	  // cout<<"sum is \n"<<ML_output_e+real_e<<endl;
	  //	  cout<<"D_star*(ML_output_e+real_e)=\n"<<D_star*(ML_output_e+real_e)<<endl;
	  //	  cout<<"D_star*(ML_output_e+real_e+Logical)=\n"<<D_star*(ML_output_e+real_e+L.transpose())<<endl;
	  //cout<<"zeromat is \n"<<one_x_n_zeromat<<endl;
	  if (D_star*(ML_output_e+real_e)==tempzeromat)
	    {
	      suc_count++;
	      if (debug==1){ cout<<"ML suc"<<endl;}
	    }
	}
      else {cout<<"sth wrong with the ML decoder, H*output_e!=syndrome"<<endl;}
	
    

      //after tansformation ML decoder:
      GF2mat syndrome_tilde(H_tilde.rows(),1);
      GF2mat zero_vec_for_syndrome(H_tilde.rows()-H.rows(),1);
      syndrome_tilde=merge_mat_vert(syndrome,zero_vec_for_syndrome);
      
      GF2mat algebraic_e_tilde(H_tilde.cols(),1);
      GF2mat	ML_output_e_tilde(H_tilde.cols(),1);
      GF2mat wrong_e(H_tilde.cols(),1);
      //  cout<<8<<endl;
      algebraic_decoder(H_tilde,syndrome_tilde,algebraic_e_tilde);
           if (debug==1)
	     {
	       cout<<"after trans algebraic_e: \n"<<algebraic_e_tilde<<endl;
	       cout<<"after trans algebraic_e+logical operator: \n"<<algebraic_e_tilde+H_tilde_star.transpose()<<endl;
	     }
  
      
      //  cout<<9<<endl;
      double exp_energy_ratio_trans=exp(energy (true,empty_tempmat, algebraic_e_tilde.transpose(), K_tilde))/(exp(energy (true,empty_tempmat, algebraic_e_tilde.transpose(), K_tilde))+exp(energy (true,empty_tempmat, algebraic_e_tilde.transpose()+H_tilde_star, K_tilde)));
      //   double exp_e_plus_c_energy=0;
      //  cout<<10<<endl;
      if (exp_energy_ratio_trans>0.5){ML_output_e_tilde=algebraic_e_tilde;wrong_e=algebraic_e_tilde+H_tilde_star.transpose();} else{ML_output_e_tilde=algebraic_e_tilde+H_tilde_star.transpose();wrong_e=algebraic_e_tilde;}
      //x cout<<"after trans ML_output_e_tilde=\n"<<ML_output_e_tilde<<endl;
      //  cout<<11<<endl;
      GF2mat trans_ML_output_e(H.cols(),1);
      GF2mat trans_ML_wrong_output_e(H.cols(),1);
      
      inverse_trans(trans_ML_output_e,ML_output_e_tilde,A);
      inverse_trans(trans_ML_wrong_output_e,wrong_e,A);
           if (debug==1)
	     {
	       cout<<"after inverse trans :\n"<<trans_ML_output_e<<endl;
	       cout<<"after inverse trans wrong e:\n"<<trans_ML_wrong_output_e<<endl;
	     }
      //cout<<12<<endl;
      //  cout<<H<<endl;
      //  cout<<trans_ML_output_e<<endl;
      if (H*trans_ML_output_e==syndrome)
	{
	  // cout<<13<<endl;
	  if (D_star*(trans_ML_output_e+real_e)==tempzeromat){after_trans_suc_count++;}
	  // else  {cout<<"(D_star*(trans_ML_wrong_output_e+real_e)=\n"<<D_star*(trans_ML_wrong_output_e+real_e)<<endl;}
	  // cout<<14<<endl;
	}
      else {cout<<"sth wrong with the ML decoder, H*output_e!=syndrome"<<endl;}
      
    }

  after_trans_theoric_suc_rate=1.0*after_trans_theoric_suc_count/max_num_cws;
  after_trans_suc_rate=1.0*after_trans_suc_count/max_num_cws;
  suc_rate=1.0*suc_count/max_num_cws;
  min_wt_suc_rate=1.0*min_wt_suc_count/max_num_cws;
}

void get_logical(const GF2mat& H, const GF2mat& D,const GF2mat& H_star, GF2mat& L)
{
  vector<GF2mat> Lo;
  GF2mat T;
  ivec perm;
  GF2mat U;
  int prerank=D.T_fact(T,U,perm);
  int temprank=0;
  int first_lo=0;
  GF2mat tempmat2=D;

  for (int i=0;i<H_star.rows();i++)
    {
      GF2mat tempmat1=H_star.get_row(i).transpose();
      tempmat2=merge_mat_vert(tempmat2,tempmat1);
      temprank=tempmat2.T_fact(T,U,perm);
      if (temprank==prerank+1){Lo.push_back(tempmat1);prerank=temprank;}
      else if(temprank!=prerank) {cout<<"sth wrong with the get_locial operator"<<endl;return;}
    }
  L=Lo[0];

  for (int j=1;j<Lo.size();j++)
    {
      L=merge_mat_vert(L,Lo[j]);
    }
  
      

  
}

//it is better to use L instead of H_star, but I havenot wrotten the function to calculate L yet.
void ML_decoder(const GF2mat& input_e, const GF2mat& D,const GF2mat& L,vector<double>& K,GF2mat& output_e){
  
  if(D.rows()>10){cout<<" D>rows()>10, it is better not to use ML decoder"<<endl;return;}
  
  double unnormalized_p=0;
 
  
  GF2mat input_e_T=input_e.transpose();
  GF2mat temp_output_e_T=input_e_T;
  double max_unnormalized_p=exp(energy (true,D, input_e_T, K));

  //for input_e:
  //start with i=1 because alpha=0 is already added
  for (int i=1;i<pow(2,D.rows()-1);i++)  
    {
      int tempi=i;
      GF2mat alpha(1,D.rows());

      // get an alpha
      for (int l=0;l<D.rows();l++)
	{
	  int tempi_red=tempi%2;
	  tempi=tempi/2;
	  alpha.set(0,l,tempi_red);
	}
      
      GF2mat tempmat1= input_e_T+alpha*D;

      max_unnormalized_p=max_unnormalized_p+exp(energy (true,D, tempmat1, K));
    }

  //calculate sum of p(e+alpha*D) over alpha for all possible e
  //go over possible e
  for (int j=0;j<L.rows();j++)
    {
      //  cout<<"j="<<j<<endl;
      //  cout<<"rows="<<H_star.rows()<<endl;
      GF2mat tempmat=L.get_row(j);
      tempmat=tempmat.transpose();

      GF2mat temp_e_T=input_e_T+tempmat;
      unnormalized_p=exp(energy (true,D, temp_e_T, K));
   
      //go over alpha
      for (int i=1;i<pow(2,D.rows()-1);i++)  
	{
	  //  cout<<"i="<<i<<endl;
	  int tempi=i;
	  GF2mat alpha(1,D.rows());

	  // get an alpha
	  for (int l=0;l<D.rows();l++)
	    {
	      int tempi_red=tempi%2;
	      tempi=tempi/2;
	      //  cout<<62<<endl;
	      alpha.set(0,l,tempi_red);
	    }
	  //	  cout<<63<<endl;
	  // cout<<temp_e_T<<endl;
	  //  cout<<D<<endl;
	  //  cout<<alpha<<endl;
	  GF2mat tempmat1= temp_e_T+alpha*D;
	  //	  cout<<64<<endl;
	  unnormalized_p=unnormalized_p+exp(energy (true,D, tempmat1, K));
	  //	  cout<<65<<endl;
	}
      if (unnormalized_p>max_unnormalized_p){max_unnormalized_p=unnormalized_p;temp_output_e_T=temp_e_T;}
      //  cout<<66<<endl;
    }
  //  cout<<67<<endl;
  output_e=temp_output_e_T.transpose();
  //  cout<<68<<endl;

}
  
//GF2mat Logical_O()

//give a output_e satisfy H*output_e=s
void algebraic_decoder(const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e) //r is the rank of H
{ 
    //get the first permutation that abs_LLR is in descending order
  int n=H.cols();
  int r=H.rows();

  GF2mat H1=H;
   
  GF2mat T;
  ivec perm2;
  GF2mat U;
  // now we have T*H1*perm2*perm2_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*e: U*e1=s1

  int rankH=H1.T_fact(T,U,perm2);
  GF2mat perm2_mat=col_permutation_matrix(perm2);
  GF2mat perm2_mat_inv=perm2_mat.inverse();
  
  // cout<<"H is \n"<<H<<endl;
  // cout<<"syndrome is \n"<<syndrome<<endl;
  
  GF2mat syndrome1=T*syndrome;
  // cout<<"syndrome1 is \n"<<syndrome1<<endl;
  
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
  GF2mat H2=U.get_submatrix(0,0,rankH-1,n-1);

  // cout<<"H2 is\n"<<H2<<endl;
  GF2mat syndrome2=syndrome1.get_submatrix(0,0,rankH-1,0);
  // cout<<"syndrome2 is\n"<<syndrome2<<endl;
  
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  // cout<<51<<endl;
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
  // cout<<"H_S is \n"<<H_S<<endl;
  GF2mat H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);
  // cout<<52<<endl;
  GF2mat HS_inv=H_S.inverse();
  // cout<<"HS_inv is\n"<<HS_inv<<endl;
  
  GF2mat e_S=HS_inv*syndrome2;
  //  cout<<"e_S is\n"<<e_S<<endl;
  GF2mat e_T(n-rankH,1);


  // cout<<53<<endl;
  for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
  //  cout<<54<<endl;
  for (int i=rankH;i<n;i++){output_e.set(i,0,e_T(i-rankH,0));}
  // cout<<55<<endl;
  output_e=perm2_mat*output_e;
  //cout<<"output_e is \n"<<output_e<<endl;
   	  	  
  if (H*output_e==syndrome){}
  else
    {
      cout<<"error H*output_e!=syndrom"<<endl;
    }  
}			 
	

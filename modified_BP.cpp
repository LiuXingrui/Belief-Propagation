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
 
  for (int i=0;i<num_of_row_reduction;i++)
    {
      int min_wt=H.cols();
      int this_row=0;
      
      for (int j=0;j<input_D.rows();j++)
	{
	  
	  bvec temp_vec=input_D.get_row(j);
	  int temp_wt=Weight(temp_vec);
	  if (temp_wt<min_wt){min_wt=temp_wt;this_row=j;}
	  
	}
      if (debug==1)
	{
	  cout<<i+1<<"th row reduction, delete "<<this_row<<"th row:"<<endl;
	}
      //perform row reduction for the column with minimum weight of D
      row_reduction(this_row,input_H, input_D, output_H,output_D,input_K,output_K,A,debug);//I havn't delete the rows and columns
 
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
	  if (i!=D.rows()-1) { output_D=output_D.get_submatrix(0,1,output_D.rows()-1,output_D.cols()-1);	}
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
	}

      
      input_H=output_H;
      input_D=output_D;
      input_K=output_K;
    }

  H_tilde=output_H;
  D_tilde=output_D;
  K_tilde=output_K;

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

  Ai.push_back(w);
  A.push_back(Ai);
  if (w>1)
    {
      Add_cols(H,D,Ai,w,debug);// now H -> H*A^T^(-1), D-> DA
      
      if (w>2)
	{
	  // cout<<"before add cols and rows"<<endl;
	  Add_cols_and_rows(H,D,Ai,w,H2,D2,debug,b);  // now H-> \frac{H*A^T^(-1)}{F}, D-> DA | DAB
	  // cout<<"after add cols and rows"<<endl;
	}
      else{H2=H;D2=D;}
      K_trans(input_K,output_K,Ai,w,b);
      // cout<<"after K_trans"<<endl;
    }
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
      int temp=1;
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
   int temp=0;
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


bool new_decoder2(GF2mat& H,GF2mat& G,GF2mat &H_tilde, GF2mat &H_tilde_star,const vector<vector<int>>& A,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,const int&  num_row_red,double& num_iter, int lmax,int& max_fail,int&syn_fail,int debug, vec &LR,int rankH,int options)
{
  int v=H.cols();
  int c=H.rows();
  int vt=H_tilde.cols();
  int ct=H_tilde.rows();

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
  GF2mat zero_vec_for_syndrome(ct-c,1);
  
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
  GF2mat syndrome_tilde(ct,1);
  syndrome_tilde=merge_mat_vert(syndrome,zero_vec_for_syndrome);
  
  mat mcv(ct,vt);   // the messages from check nodes to variable nodes, which are p0/p1
  mat mvc(ct,vt);   // the messages from variable nodes to check nodes
  mcv.zeros();
  mvc.zeros();
    
  initialize_massages( mcv,mvc, H_tilde); //initialize to all-1 matrix

  GF2mat output_e(v,1);
  GF2mat output_e_tilde(vt,1);
      
  for (int l=1;l<=lmax;l++)
    {
      //  cout<<l<<endl;
      quan_s_update(checks,errors, mcv,mvc,syndrome_tilde,pv_dec, ct, vt,output_e_tilde,LR,1);	
	    
      if (H_tilde*output_e_tilde==syndrome_tilde)
	{
	  //inverse_trans(output_e,output_e_tilde,A);
	  cout<<"ok"<<endl;
	  return true;
	  if(H*(output_e+real_e)==zero_rvec2)
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

#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <chrono>
#include <thread>

#ifndef BP
#define BP
#include"BP.h"
#endif


using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

void initialize_checks (const GF2mat &H, nodes checks[], int & E){
  int r=H.rows();
  int c=H.cols();
  
  for (int i=0; i<r;i++)
    {
    checks[i].degree=0;
    for (int j=0; j<c;j++)
      { 
      if (H(i,j)==1)
	{
	checks[i].degree++;
	(checks[i].neighbors).ins(0,j);
	E++;
      }     
    }    
  }
}


//find the neighbors of variable(error) nodes
void initialize_errors(const GF2mat &H, nodes errors[]){

    int r=H.rows();
    int c=H.cols();
    int index=0;
 
    for (int i=0; i<c;i++)
      {
      errors[i].degree=0;
      for (int j=0; j<r ;j++)
	{
	if (H(j,i)==1)
	  {
	  errors[i].degree++;
	  (errors[i].neighbors).ins(0,j);
	  index++;	
	}	
      }
    }
}

//initialize massages to all-one matrices
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H){
    int r=H.rows();
    int c=H.cols();
    for (int i=0;i<r;i++)
      {
      for (int j=0;j<c;j++)
	{
	    if (H(i,j)==1){
	      mcv.set(i,j,1);
	      mvc.set(i,j,1);	 
	     }
	}	     
    }  
}

//parallel schedule for classical BP
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e){
  double ipr=(1-p)/p;
   
  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        

	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr);	
      }
    }
    //update all c-to-vj massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr=ipr;

   for (int i=0;i<vj_degree;i++)
     {
       int cnode=(errors[j].neighbors)(i);
       update_ci_to_vj( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));

       final_pr=final_pr*mcv(cnode,j);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;
      output_e.set(j,0,final_pr<1? 1:0);        
    }  
}


//the sequential schedule 
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat &mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e)
{
    double ipr=(1-p)/p; //initial p0/p1
    double  final_pr;  //final p0/p1
    
    //fix j,  for every v_j, do the following:
    for (int j=0;j<v;j++)
      {
	 int vj_degree=errors[j].degree;
	 final_pr=ipr;

	 //ci is the ith neibor of vj:
	 for (int i=0;i<vj_degree;i++)
	   {
	     //  update the  v_j to c_i massage:      
	     int ci=(errors[j].neighbors)(i);
	     update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);          
	   }
	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));	   // update all  c_i to v_j massages         
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 output_e.set(j,0,final_pr<1? 1:0);      
      }  
}

//update vi to cj massage:
inline void update_ci_to_vj(const nodes checks[],const nodes errors[],mat& mcv,mat& mvc,int i,int j,bin s)
{
   int c_degree=checks[i].degree;
   double temp=1.0;

   for (int kk=0;kk<c_degree;kk++)
     {
      int vk=(checks[i].neighbors)(kk);
      if (vk!=j)
	{
	  temp=temp*(mvc(i,vk)-1)/(mvc(i,vk)+1);	  
	}
     }

     mcv.set(i,j,s==0? (1+temp)/(1-temp):(1-temp)/(1+temp));
     
     if (isinf(mcv(i,j)))
       {
	 //cout<<"inf"<<endl;
	 mcv.set(i,j,1);
       }
     
}

//update vj to ci massage:
inline void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr,double alpha)
{ 
   int v_degree=errors[j].degree;
   mvc.set(i,j,ipr);

   for (int kk=0;kk<v_degree;kk++)
     {
      int ck=(errors[j].neighbors)(kk);
      if (ck!=i)
	{
	  mvc.set(i,j,mvc(i,j)*pow(mcv(ck,j),1.0/alpha));
	}
      else
	{
	   mvc.set(i,j,mvc(i,j)*pow(mcv(ck,j),1.0-1.0/alpha));
	}
   
     }   
}

//clasical decoding
int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,vec& pv,int debug)
{
  int n=v; 
  //if no error, break
  double p=pv[0];
  int wt=0;
  GF2mat real_e(v,1);    
  wt=error_channel(real_e, pv);
  //GF2mat zero_vec(v,1);
  if (wt==0)
    {
      return 1;
    }
  
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{	    
	  er=er+ distance(zero_mat2, real_e, n); // for calculating bit error rate after decoding
	  return 0;    //failed		   
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
	  if (debug==0)
	    {
	      p_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	    }
	  else
	    {
	      s_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	    }
	  if (H*output_e==syndrome)
	    {
	      num_iter=num_iter+l;
	     
	      if(output_e==real_e)
		{
		  //cout<<"success! iteration number="<<l<<endl;
		  return 1;
		  
		  // num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	      else
		{
		  //cout<<"error!,get the wrong e"<<endl;
		  er=er+ distance(output_e, real_e, n);	  
		  return 0;		       
		}	    	  
	    }
	  if(l==lmax)
	    {
		er=er+ distance(output_e, real_e, n);
		cout<<"real_e"<<endl;
        
	      return 0;	   	 
	    }
	}
 }

  //the distance between 2 cws
int distance(const GF2mat &output_e, const GF2mat &real_e,int n)
{ 
  int er=0;
  for (int i=0;i<n;i++)
    {
    if( (output_e(i,0)!=real_e(i,0)))
      {
      er++;
      }
    }
  return er;
}


//read a parity check matrix
GF2mat read_matrix (int& n,int &r, string & file_name)
{
  ifstream parity_check(file_name);
  string line;
  int row_ind=0;
  int temp;
  getline(parity_check, line);
  istringstream iss(line);
  
  if(  iss>>n>>r) {}
  else
    {
      cout<<"did not open the right parity check file"<<endl;
      
    }
 
  GF2mat H(r,n);
  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp){
	if (temp>=1&&temp<=n&&row_ind<r)
	  {
	    H.set(row_ind,temp-1,1);
	    // cout<< H(row_ind,temp-1)<<endl;
	  }
	else
	  {
	    cout<<"the format of the parity check is wrong, the first element is 1 rathar than 0"<<endl;	    
	  }	
      }
      row_ind++;
    }
  parity_check.close();
  // cout<<H<<endl;
  return H;  
}


void write_matrix(string file_name, GF2mat &H)
{
  ofstream Hx;
  Hx.open (file_name,ios::trunc);
  int n=H.cols();
  int r=H.rows();

  Hx<<n<<" "<<r<<endl;

  for (int j=0;j<r;j++)
    {
      for (int i=0; i<n;i++)
	{
	  if (H(j,i)!=0)
	    {
	      Hx<<i+1<<" ";
	    }	  
	}
      Hx<<endl;
    }
  Hx.close();
}


GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right)
{
  if (left.rows()!=right.rows())
    {
      cout<<"the 2 matrices cannot merge horizontally"<<endl;
      GF2mat error(1,1); 
      return error;
    }

  else
    {
      int r=left.rows();
      int c=left.cols()+right.cols();
      int c1=left.cols();
      int c2=right.cols();
      GF2mat m(r,c);

      for (int i=0;i<r;i++)
	{
	  for (int j1=0;j1<c1;j1++){
	    m.set(i,j1,left(i,j1));
	  }
	  for (int j2=0;j2<c2;j2++){
	    m.set(i,c1+j2,right(i,j2));
	  }	  
	}
       return m;
    }
}

//the bsc error channel
int error_channel(GF2mat &cw, const vec &p)
{
  double temp2;
  bin one=1;
  int temp=0;
  int r=cw.rows();
  if (cw.rows()!=p.size())
    {
      cout<<"the size of p and cw do not match"<<endl;
    }
  else
    {
    for (int i=0;i<r;i++)
      {
	temp2=randu();
	if(temp2<p[i])
	  {	 
	    cw.set(i,0,cw(i,0)+one);
	    temp++;
	  }	
      }
    }
  return temp;
}

// error channel for fixed weight errors
void error_channel2(GF2mat &error, int wt)
{ 
  double temp2;
  bin one=1;
  int n=error.rows();

    for (int i=0;i<wt;i++)
      {    
	temp2=randi(0,n-1);
	error.set(temp2,0,one);	      
      }
    if (s_weight(error)!=wt)
      
      { GF2mat error2(n,1);
	error=error2;
	//cout<<"!=wt, try again"<<endl;
	error_channel2(error,wt);
      }
}

//weight of a row vector
int weight(const GF2mat &cw)
{
  int n=cw.cols();
  int wt=0;
  for (int i=0;i<n;i++)
    {
      if(cw(0,i)==1){wt++;}
    }
  return wt;
}

//weight of a column vector
int s_weight(const GF2mat &s)
{
  int k=s.rows();
  int wt=0;
  for (int i=0;i<k;i++)
    {
      if(s(i,0)==1){wt++;}
    }
  return wt;
}


//get the error rate distribution
void pro_dist(double pmin,double pmax, vec& pv)
{ 
  double pdiff=pmax-pmin;
  int pvsize=pv.size();
  double temp;
 
  for (int i=0;i<pvsize;i++)
    {
      temp=pdiff*randu();
      pv(i)=pmin+temp;
    }
}


bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,int& max_fail, int&syn_fail,  int debug, vec &LR,int rankH,int &OSD_suc,double alpha,double lambda)
{
  int v=H.cols();
  int c=H.rows();

  int wt_real_e=0;
  GF2mat real_e(v,1); 

  if (wt==0)
    {
      wt_real_e=error_channel(real_e, pv);
      if (wt_real_e==0)
	{	  
	  return true;
	}
    }
  else
    {
      error_channel2(real_e,wt);
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
	  if ((debug/2)%2==1)
	    {
	      quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,alpha);
	    }
	  else
	    {
	      quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,alpha);
	    }
	    
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
      	            
      if (debug%2==1)
	{	  
	  mcv.zeros();
	  mvc.zeros();
	  initialize_massages( mcv,mvc, H);
	  vec pv2(v);
	  pro_dist(pmin,pmax,pv2);
	  GF2mat output_e2(v,1);
	  GF2mat syndrome2=H*(real_e+output_e);
	  for (int l=1;l<=lmax;l++)
	    {	   
	      if ((debug/2)%2==1)
		{
		  quan_p_update(checks,errors, mcv,mvc,syndrome2,pv2, c, v,output_e2,LR,alpha);
		}
	      else
		{
		  quan_s_update(checks,errors, mcv,mvc,syndrome2,pv2, c, v,output_e2,LR,alpha);
		}
	      if (H*output_e2==syndrome2)
		{		  
		  if(G*(output_e+real_e+output_e2)==zero_rvec2)
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
	     output_e=output_e+output_e2;
	}
   
      if ((debug/8)%2==1)
	{
	  OSD(LR,H,syndrome,output_e);
	  if(G*(output_e+real_e)==zero_rvec2)
	    {
	      // cout<<"OSD suc2.1"<<endl;
	      OSD_suc++;
	      return true;
	    }	 
	}
       if ((debug/4)%2==1)
	 {
	   cout<<"output e is \n"<<endl;
	   err_pos2(output_e);
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	 }
      max_fail++;
      return false;
 }

//ordered statistical decoder
void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e) //r is the rank of H
{ 
    //get the first permutation that abs_LLR is in descending order
  int n=LR.length();
  int r=H.rows();
  vec LLR(n); //the sort function gives a ascending  order, LR=p0/p1
  for (int i=0;i<n;i++)
    {
      if (LR(i)>=0) 
	{
	  LLR(i)=log(LR(i));
	}
      else
	{
	  cout<<"negative LR!"<<endl;
	  cout<<LR(i)<<endl;
	}
    }
  
  //now we have H*perm1*perm1_inv*e=s, and def:H1=H*perm1
  ivec perm1=sort_index(LLR);
  //for (int i=0;i<n;i++){perm1(i)=i;}
  GF2mat H1=H;
  H1.permute_cols(perm1,0);
  GF2mat perm1_mat=col_permutation_matrix(perm1);
  GF2mat perm1_mat_inv=perm1_mat.inverse();
   
  GF2mat T;
  ivec perm2;
  GF2mat U;
  // now we have T*H1*perm2*perm2_inv*perm1_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*perm1_inv*e:

  int rankH=H1.T_fact(T,U,perm2);
  GF2mat perm2_mat=col_permutation_matrix(perm2);
  GF2mat perm2_mat_inv=perm2_mat.inverse();
      
  GF2mat syndrome1=T*syndrome;
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
  GF2mat H2=U.get_submatrix(0,0,rankH-1,n-1);
  GF2mat syndrome2=syndrome1.get_submatrix(0,0,rankH-1,0);
  
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
  GF2mat H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);
  GF2mat HS_inv=H_S.inverse();
  
  GF2mat e_S=HS_inv*syndrome2;  
  GF2mat e_T(n-rankH,1);
  GF2mat new_e_S;
  int wt=s_weight(e_S)+1;
  //cout<<"original wt "<<wt<<endl;
  int temp1;
  int temp2=-1;
  int temp3=-1;
  int lambda=min(60,n-rankH);
  
  for (int i=0;i<n-rankH;i++)
    {
      //cout<<i<<": wt is "<<wt<<endl;
      GF2mat new_e_T=e_T;
      new_e_T.set(i,0,1);
      new_e_S=HS_inv*e_S+HS_inv*H_T*e_T;
      temp1=s_weight(new_e_S);
      if (temp1+1<wt)
	{
	  wt=temp1+1;
	  temp2=i;
	}
    }

  for (int i=0;i<lambda;i++)
    {
      for (int j=0;j<lambda;j++)
	{
	  GF2mat new_e_T=e_T;
	  new_e_T.set(i,0,1);
	  new_e_T.set(j,0,1);		  
	  new_e_S=HS_inv*e_S+HS_inv*H_T*e_T;
	  temp1=s_weight(new_e_S);

	  if (temp1+2<wt)
	    {
	      // cout<<"another e_S"<<endl;
	      wt=temp1+2;
	      // cout<<"new wt is "<<wt<<endl;
	      temp2=i;
	      temp3=j;
	    }
	}
    }
 
      if (temp2!=-1&&temp3==-1)
	{
	  e_T.set(0,temp2,1);
	  e_S=e_S+HS_inv*H_T*e_T;
	}
      else if (temp2!=-1&&temp3!=-1)
	{
	  e_T.set(0,temp2,1);
	  e_T.set(0,temp3,1);
	  e_S=e_S+HS_inv*H_T*e_T;
	}
           
      for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
      for (int i=rankH;i<n;i++){output_e.set(i,0,e_T(i-rankH,0));}
      output_e=perm1_mat*perm2_mat*output_e;
   	  	  
      if (H*output_e==syndrome){}
      else
	{
	  cout<<"error H*output_e!=syndrom"<<endl;
	}  
}			 
	
//get a column permutation matrix
GF2mat col_permutation_matrix(ivec & perm)
{
  int n=perm.length();
  GF2mat p(n,n);
  for ( int i=0;i<n;i++)
    {
      p.set(perm(i),i,1);
    }
  return p;
}

/*
GF2mat col_permutation_matrix_s(ivec & perm)
{
  int n=perm.length();
  GF2mat p(n,n);
  for ( int i=0;i<n;i++)
    {
      p.set(perm(i),i,1);
    }
  return p;
}
*/

void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha)
{  
  double ipr;
  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;      
      //update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        ipr=(1-pv[vnode])/pv[vnode];
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr,alpha);	
      }
    }
    //update c-to-vj massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr=(1-pv[j])/pv[j];
   
   for (int i=0;i<vj_degree;i++)
     {
       int cnode=(errors[j].neighbors)(i);
       update_ci_to_vj( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));
       final_pr=final_pr*pow(mcv(cnode,j),1.0/alpha);
      }   
   LR(j)=final_pr;
   output_e.set(j,0,final_pr<1? 1:0);   
    }  
}


void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha)
{  
    double ipr;
    double final_pr;
   for (int j=0;j<v;j++)
      {
	 int vj_degree=errors[j].degree;
	 ipr=(1-pv[j])/pv[j];
	 final_pr=ipr;
	 //ci is the ith neibor of vj:
	 for (int i=0;i<vj_degree;i++)
	   {
	     //  update the  v_j to c_i massage:      
	     int ci=(errors[j].neighbors)(i);
	     update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);          
	   }	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));	   // update all  c_i to v_j massages         
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 LR(j)=final_pr;
	 output_e.set(j,0,final_pr<1? 1:0);      
      }  
}

GF2mat get_gen(const GF2mat &H){

  GF2mat HT=H.transpose();
  GF2mat T,U;
  ivec P;
  int Hrank= HT.T_fact(T,U,P);
  int r=H.rows();
  int n=H.cols();
  GF2mat G=T.get_submatrix(Hrank,0,n-1,n-1);
  return G;
}
/*
void dense_to_sparse(GF2mat &G,GF2mat& G_s)
{
  
  int r=G.rows();
  int c=G.cols();
  for (int i=0;i<r;i++)
    {
      for (int j=0;j<c;j++)
	{
	  if (G(i,j)==1){G_s.set(i,j,1);}	
	}
    }
}
*/
void err_pos1(const GF2mat &error)
{
  cout<<"wt="<<s_weight(error)<<endl;
    int n=error.rows();
    for (int i=0;i<n;i++)
      {
	if (error(i,0)==1)
	  {
	    cout<<i<<" ";
	  }
      }
    cout<<endl;
}
void err_pos2(const GF2mat &error)
{
  cout<<"\n";
  int n=error.rows();
  int d=sqrt(n);
  if (d*d==n)
    {
      for (int r=0;r<d;r++)
	{
	  cout<<r<<"th row: ";
	  for (int c=0;c<d;c++)	    
	    {
	      if (error(r*d+c,0)==1){cout<<1;}
	      else {cout<<".";}
	    }
	  cout<<endl;
	}
    }
    cout<<"\n";
}

int GF2mat_rank(const GF2mat& H){

  GF2mat T,U;
  ivec P;
  return H.T_fact(T,U,P);	
}

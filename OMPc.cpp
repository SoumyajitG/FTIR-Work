#include<iostream>
#include<cmath>
#include <stdlib.h>
#include<memory>
#include"mex.h"
#include "C:\Users\Lenovo\Documents\Visual Studio 2012\Projects\eigenpinv\eigenpinv\Eigen\Eigen" 
#include<windows.h>
using namespace std;
using namespace Eigen;
void ginv(double * a,double *AZheng ,int m, int n) 
{
	 MatrixXf A(m,n);   //动态矩阵，建立3行4列。  
     for (int i = 0;i<m;i++)
         for(int j = 0;j<n;j++)
             A(i,j) = a[i*n+j];
    // cout<<A<<endl;
JacobiSVD<MatrixXf> svd(A, ComputeFullU | ComputeFullV);
VectorXf sig = svd.singularValues();
//cout<<"sig"<<sig<<endl<<"LL";
double  pinvtoler=1.e-6; // choose your tolerance wisely!
     for ( long i=0; i<A.cols(); ++i) {
        if ( sig(i) > pinvtoler )
           sig(i)=1.0/sig(i);
       else sig(i)=0;
     }
	 int col = A.cols(),row =A.rows();
	 int d = min(col,row);
	 MatrixXf sigu = MatrixXf::Zero(col,row);

	 for ( long i=0; i<d; ++i) {
        sigu(i,i) = sig(i);
     }
	 MatrixXf V = svd.matrixV();
	 MatrixXf U = svd.matrixU().transpose();
	 MatrixXf ans = V*sigu;
	 MatrixXf an = ans*U;
	//cout<<sigu<<endl;
      for (int i = 0;i<n;i++)
         for(int j = 0;j<m;j++)
             AZheng[i*m+j] = an(i,j);


	
}
int absmaxpos(double *a,int M)
{  
	int pos=0;
	double max=a[0];
	for(int i=0;i<M;i++)
	{
		a[i]=fabs(a[i]);
		if(max<=a[i])
		{	
			max=a[i];
			pos=i;
			
		}
	}
	return pos;
}
void trmul(double *a,double *b,int m,int n,int k,double *c)
{ 
	int i,j,l,u;
	for (i=0; i<=m-1; i++)
		for (j=0; j<=k-1; j++)
		{ 
			u=i*k+j;
			c[u]=0.0;
			for (l=0; l<=n-1; l++)
				c[u]=c[u]+a[i*n+l]*b[l*k+j];
		}
}
int Trans(double *matrix,int M,int N)
{
	int i,j;
	int tol=M*N*sizeof(double);
	double *temp=new double[tol];
	memset(temp,0,tol);
	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
			temp[i*M+j]=matrix[j*N+i];
		}
	}
	memcpy(matrix,temp,tol);
	delete[] temp;
	return 0;
}
void print(double *a, int row,int col){//hang zhu xhu
    for(int i = 0;i<row;i++)
    {
        for (int j = 0;j<col;j++)
        {
            printf("%f\t",a[i*col+j]);
        }
        cout<<endl;
    }
}
            
void OMP(double *b,double *A,int M,int N,int K,double *result)  //T1为车辆
{//A  M*N  
    int i,j;
    int pos;
	double norm;
	for(int i = 0;i<N;i++){
		result[i] = 0;
	}
	double *r =new double[M];
	double *x_T = new double[K];
	memcpy(r,b,sizeof(double)*M);
    double *Ar=new double[N];
	int *indx_set = new int[K];
	double*A_T = new double[M*K];
	double*A_T_nonorth = new double[M*K];
	trmul(r,A,1,M,N,Ar);
 //    print(Ar,N,1);
	for (int kk=0;kk<K;kk++)
    {
     //   cout<<"round :"<<kk<<endl;for(int i = 0;i<N;i++)cout<<"AR:"<<Ar[i];cout<<endl;
		int ind_new = absmaxpos(Ar,N);
		indx_set[kk] = ind_new;//cout<<"ind_new:"<<ind_new;cout<<endl;
		//sort(indx_set,K);
		double* atom_new = new double[M];
		for(int i = 0;i<M;i++){
			atom_new[i] = A[i*N+ind_new];
			A_T_nonorth[i*K+kk] = A[i*N+ind_new];
		}
       // print(A,M,N);
    //    print(A_T_nonorth,M,K);
		for(int j = 0;j<kk;j++){
			int coff = 0;
			for(int l = 0;l<M;l++){
				coff += atom_new[l]*A_T[l*K+j];
			}
			for(int l = 0;l<M;l++){
				atom_new[l] -= coff*A_T[l*K+j];
           //     cout<<atom_new[l]<<endl;
			}
		}
		norm = 0;
		for(int i=0;i<M;i++){
			norm += atom_new[i]*atom_new[i];
		}
		norm = sqrt(norm);
		for(int i=0;i<M;i++){
			atom_new[i] /= norm;
		}
     //   print(atom_new,M,1);
		for(int i=0;i<M;i++){
			A_T[i*K+kk] = atom_new[i];
		}
      //  print(A_T,M,K);
		double* tmpA_T = new double[M*(kk+1)];
		for(int i =0 ;i<M;i++){
			for(int j = 0;j<kk+1;j++){
				tmpA_T[i*(kk+1)+j] = A_T[i*(K)+j];
			}
		}
     //   print(tmpA_T,M,kk+1);
      //  Trans(tmpA_T,M,kk+1);
		trmul(b,tmpA_T,1,M,kk+1,x_T);
    //    for(int i = 0;i<kk+1;i++)cout<<"HELLO!:"<<x_T[i]<<endl;
		for(int i = 0;i<kk+1;i++){
			result[indx_set[i]] = x_T[i];
		}
		double *tmp = new double[M];
		//trmul(tmpA_T,x_T,M,kk+1,1,tmp);
        Trans(tmpA_T,M,kk+1);
		trmul(x_T,tmpA_T,1,kk+1,M,tmp);
		
		for(int i = 0;i<M;i++){
			r[i] = b[i] - tmp[i];//cout<<r[i];
		}//cout<<endl;
		if (kk<K){
			trmul(r,A,1,M,N,Ar);
		}
        
	}
	double *pinv = new double[K*M];
   // FILE* fp = fopen("test.txt","a+");
 // for(int m = 0;m<M;m++){
  //    cout<<b[m]<<endl;
  //    for(int n=0;n<K;n++)
   //       fprintf(fp,"%f\t",A_T_nonorth[m*K+n]);
    //  fprintf(fp,"\n");
  //}
   // fclose(fp);
	ginv(A_T_nonorth,pinv,M,K);
 //   print(pinv,K,M);
    trmul(pinv,b,K,M,1,x_T);
	for(int i = 0;i<K;i++){
		result[indx_set[i]] = x_T[i];
	}
    delete []r;delete []x_T;delete []Ar;delete []indx_set;delete []A_T;delete []A_T_nonorth;
}



void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) //OMP(y,d,t)
{ 
   
	double *Y,*D,*sc,*T; 
	int M,N,K; 
	Y=mxGetPr(prhs[0]); 
	M=mxGetM(prhs[1]); 
	K=mxGetN(prhs[1]);
	N=mxGetN(prhs[0]);
	D=mxGetPr(prhs[1]);
	T = mxGetPr(prhs[2]);
	double *omptmp = (double*)malloc(sizeof(double)*K);
	plhs[0]=mxCreateDoubleMatrix(K,N,mxREAL); 
	sc=mxGetPr(plhs[0]); 
	double*y = new double[M];
	double* Dtmp = new double [M*K];
	for(int i=0;i<N;i++){
		for(int j = 0;j<M*K;j++){
			Dtmp[j] = D[j];
        }
        Trans(Dtmp,K,M);
		for(int j = 0;j<M;j++){
			y[j] = Y[i*M + j];
        }

		OMP(y,Dtmp,M,K,int(*T),omptmp);

		for (int j = 0;j<K;j++){
			sc[i*K + j] = omptmp[j];
           // cout<<omptmp[j]<<endl;
           // sc[i*K + j] =0;
		}
    }
//    cout<<"gg"<<endl;
	delete []Dtmp;
	delete []y;
    
}









#include<iostream>
#include<cmath>
#include <stdlib.h>
#include<memory>
#include"mex.h"
using namespace std;
void ginv(double * a,double*AZheng ,int m, int n) ;
int LineSimpleJuZhen(double * a,int m2,int n2);//һ�λ���
int From_aa_to_a1(double * aa,int mm,int nn,double * a1,int j1);//����aa[mm*nn]�õ���һ����һ�еľ���a1[(mm-1)*nn-(j1+1)]
//ת�Ⱦ���,��a��Ϊat,����aΪPointNum*two�� 
int AtoAT(double a[],double at[],int PointNum, int two);//רָ����,aΪ PointNum*two�� 
//����˷���C=A*B��AΪPointNum*two,BΪtwo*two_B��CΪPointNumb*two_B
void ChengFa(double *const A, double *const B, double *const C,int PointNum,int two,int two_B);
//��������
void ContraryMatrix(double *const pMatrix, double *const _pMatrix, const int &dim);
void ginv(double * a,double*AZHENG, int m, int n) 
{
	double * b=new double[m*n];///��b������a��ԭʼ���� ������a����д�� ԭʼ���ݶ�ʧ
	for(int i100=0;i100<m*n;i100++)
	{
		b[i100]=a[i100];
	}
	///(һ����Ϊ���ݾ���aa


	int j1;//��һ�εķ��������ڵ�����
	int mm;//����ѭ���е�aa�����Ԫ�ظ�������
	int nn;//����ѭ���е�aa�����Ԫ�ظ�������
	int Count=0;//������Count�μ򻯴���
	int JuZhenNotZero=0;//�ж��Ƿ�Ϊ�����,0��ʾ�����������
	mm=m;
	nn=n;
	double * aa=new double[mm*nn];//����aaʵ��ÿ��ѭ������������Ĵ�С��һ�£���Ϊaa��Ԫ�ش�С�ɱ�
	
	
	aa=a;

	j1=LineSimpleJuZhen(aa,mm,nn);//��aa[mm][nn]���е�һ�λ���,j1Ϊ����ֵ������ ��aa Ϊ�����ľ���
	Count=Count+1;

	//��ʣ�µ�A1 
	double * aaa1=new double[(mm-Count)*(nn-(j1+1))];
	JuZhenNotZero = From_aa_to_a1(aa,mm,nn,aaa1,j1);

	

	double * ptaaa1;//���ڴ���aaa1�ĵ�ַ��aa1
	ptaaa1=aaa1;
	int j2=0;
	while (JuZhenNotZero==1)///aaa1!=0����
	{
		mm=m-Count;
		nn=n-(j1+Count);
		double * aa1=new double[mm*nn];
		aa1=ptaaa1;
		j2=LineSimpleJuZhen(aa1,mm,nn);
		
  
	    //ԭʼ��aa����и���

		for(int i5=Count;i5<m;i5++)
			for(int j5=j1+Count;j5<n;j5++)
			{
				aa[i5*n+j5]=aa1[(i5-Count)*nn+(j5-(j1+Count))];
			}

		Count=Count+1;//��aa���е�Count�α任����
		if((m-Count)*(nn-(j2+1))==0)
		{
			JuZhenNotZero=0;//�����Ѿ���n*1��1*n�ˣ���ʾ�Ѿ���������ѭ����
		}
		else
		{
	       double * aaa1=new double[(m-Count)*(nn-(j2+1))];//���Ӷ�(m-Count)*(nn-(j2+1))�Ƿ�Ϊ0���ж�
	       JuZhenNotZero=From_aa_to_a1(aa1,mm,nn,aaa1,j2);
		   ptaaa1=aaa1;



	       j1=j1+j2;
		  
		}
	}



	double kk=0;
	int j6=0;
	int jjj;//��jjj��Ϊ��һ��������
	for(int i6=0;i6<Count;i6++)
	{
		j6=0;
		do
		{
			kk=aa[i6*n+j6];
			j6=j6+1;
		}
		while(kk==0);//�ҵ�ÿһ�е�һ��������
		jjj=j6-1;


		for(j6=0;j6<n;j6++)
		{
			aa[i6*n+j6]=aa[i6*n+j6]*1.0/kk;

		}
		//����������ĸ��л�Ϊ0




    }

	double kkk=0;
	int jjjj=0;
	int j66;
	double kkkk=0;

	for(int i8=1;i8<Count;i8++)
	{
		j66=0;
		do
		{
			kkk=aa[i8*n+j66];
			j66=j66+1;
		}
		while(kkk==0);//�ҵ�ÿһ�е�һ��������
		jjjj=j66-1;
		for(int i9=0;i9<i8;i9++)//���л�
		{
			kkkk=aa[i9*n+jjjj]/1.0;
			for(int j9=0;j9<n;j9++)
			{
				

				aa[i9*n+j9]=aa[i9*n+j9]-aa[i8*n+j9]*kkkk;
			}
		}



	}

	double * G=new double [Count*n];
	double * F=new double [m*Count];
	int * ArrayCount=new int [Count];//���aa�еڼ����Ǵ��о�һ��1
	double kA;
	int j67;
	int jjjjj;

	for(int i11=0;i11<Count;i11++)//�õ�ArrayCount����
	{
		j67=0;
		do
		{
			kA=aa[i11*n+j67];
			j67=j67+1;
		}
		while(kA==0);//�ҵ�ÿһ�е�һ��������
		jjjjj=j67-1;

		ArrayCount[i11]=jjjjj;


	}

	for(int i10=0;i10<Count;i10++)
		for(int j10=0;j10<n;j10++)
			G[i10*n+j10]=aa[i10*n+j10];
//	DisplayMatrix(G,Count,n);
	//�õ�F//�˴�����aa�������Ѿ����ı��ˣ���b�����ݶ�F����
	for(int i13=0;i13<m;i13++)
		for(int j13=0;j13<Count;j13++)
		{
			F[i13*Count+j13]=b[i13*n+ArrayCount[j13]];

		}

	////	��A�Ĺ��������ʼ 

	double * GT=new double [n*Count];
	double * FT=new double [Count*m];

	AtoAT(G,GT, Count,  n);
	AtoAT(F,FT,m, Count);


	double * Cheng_G_GT=new double [Count*Count];
	double * Cheng_FT_F=new double [Count*Count];

	ChengFa(G,GT,Cheng_G_GT,Count,n,Count);
	ChengFa(FT,F,Cheng_FT_F,Count,m,Count);

	double * ni_Cheng_G_GT=new double [Count*Count];
	double * ni_Cheng_FT_F=new double [Count*Count];

    ContraryMatrix(Cheng_G_GT,ni_Cheng_G_GT,Count);
	ContraryMatrix(Cheng_FT_F,ni_Cheng_FT_F,Count);

	double * Cheng_GT_niG=new double [n*Count];
	ChengFa(GT,ni_Cheng_G_GT,Cheng_GT_niG,n,Count,Count);

	double * Cheng_GT_niG_niF=new double [n*Count];

	ChengFa(Cheng_GT_niG,ni_Cheng_FT_F,Cheng_GT_niG_niF,n,Count,Count);

	ChengFa(Cheng_GT_niG_niF,FT,AZHENG,n,Count,m);
}



int  LineSimpleJuZhen(double * a,int m2,int n2)//һ�λ���//����ֵΪ��һ�������е�����
{
	int j=0;
	int i;
	int break1=0;
	

	do//d�õ���һ�� ������a��i]��j]
	{
		
		for(i=0;i<m2;i++)
		{
			if(a[i*n2+j]!=0)
			{
				break1=1;
				break;
			}
		}
		if(break1==1)
			break;
		j=j+1;
	}
	while (j<n2);

	if(i!=0)
	{
		for(int j1=0;j1<n2;j1++)
		{
			double pp1;
			pp1=a[0*n2+j1];
			a[0*n2+j1]=a[i*n2+j1];
			a[i*n2+j1]=pp1;
		}
	}

	double k;
	for(int i4=1;i4<m2;i4++)
	{
		k=(0-a[i4*n2+j])/1.0/a[0*n2+j];//���ȵõ�ϵ���������ı�a��i4����j]��ֵ
		for(int j4=0;j4<n2;j4++)
		{
			a[i4*n2+j4]=a[i4*n2+j4]+k*a[0*n2+j4];
		}
	}

	return j;

}


int From_aa_to_a1(double * aa,int mm,int nn,double * a1,int j1)//����aa[mm*nn]�õ���һ����һ�еľ���a1[(mm-1)*nn-(j1+1)]
{
	int JuZhenNotZero=0;
	for(int i2=0;i2<mm-1;i2++)
		for(int j2=0;j2<(nn-(j1+1));j2++)
		{

			a1[i2*(nn-(j1+1))+j2]=aa[(i2+1)*nn+(j2+(j1+1))];//��a����A1
			//�ж��Ƿ�Ϊ�����
			if(a1[i2*(nn-(j1+1))+j2]!=0)//ֻҪ��һ�� ��Ϊ�㼴��
			{
				JuZhenNotZero=1;//1��ʾ�����Ƿ������
			}
		}

	return JuZhenNotZero;

}

//ת�Ⱦ���,��a��Ϊat,����aΪPointNum*two�� 
int AtoAT(double a[],double at[],int PointNum, int two)//רָ����,aΪ PointNum*two�� 
{
	for(int i=0;i<two;i++)
		for(int j=0;j<PointNum;j++)
		{
			at[i*PointNum+j]=a[j*two+i];
		}
	return 1;
}


//����˷���C=A*B��AΪPointNum*two,BΪtwo*two_B��CΪPointNumb*two_B
void ChengFa(double *const A, double *const B, double *const C,int PointNum,int two,int two_B)
{
	double sub=0.0;
     for (int i=0; i<PointNum; i++)
     {
         for (int j=0; j<two_B; j++)
         {
		   sub=0.0;
           for(int k=0;k<two;k++)
		   {
			   sub=sub+A[i*two+k]*B[k*two_B+j];
		   }
		   C[i*two_B+j]=sub;
         } 
     }
}



//��pMatrix������󣬲������ھ���_pMatrix��
void ContraryMatrix(double *const pMatrix, double *const _pMatrix, const int &dim)
{ 
     double *tMatrix = new double[2*dim*dim];
     for (int i=0; i<dim; i++)
	 {
         for (int j=0; j<dim; j++)
             tMatrix[i*dim*2+j] = pMatrix[i*dim+j];        
     }

     for (int i=0; i<dim; i++)
	 {
         for (int j=dim; j<dim*2; j++)
             tMatrix[i*dim*2+j] = 0.0;

         tMatrix[i*dim*2+dim+i]   = 1.0;        
     }
     //Initialization over!
     for (int i=0; i<dim; i++)//Process Cols
     {
         double base = tMatrix[i*dim*2+i];

         for (int j=0; j<dim; j++)//row
         {
             if (j == i) continue;
             double times = tMatrix[j*dim*2+i]/base;
             for (int k=0; k<dim*2; k++)//col
             {        
                 tMatrix[j*dim*2+k] = tMatrix[j*dim*2+k] - times*tMatrix[i*dim*2+k];
             }
         }        

         for (int k=0; k<dim*2; k++)
		 {
             tMatrix[i*dim*2+k] /= base;
         }
     }
     for (int i=0; i<dim; i++)
     {
         for (int j=0; j<dim; j++)
             _pMatrix[i*dim+j] = tMatrix[i*dim*2+j+dim];        
     }    
     delete[] tMatrix;
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

void OMP(double *b,double *A,int M,int N,int K,double *result)  //T1Ϊ����
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
	for (int kk=0;kk<K;kk++)
    {
		int ind_new = absmaxpos(Ar,N);
		indx_set[kk] = ind_new;
		//sort(indx_set,K);
		double* atom_new = new double[M];
		for(int i = 0;i<M;i++){
			atom_new[i] = A[i*N+ind_new];
			A_T_nonorth[i*M+kk] = A[i*N+ind_new];
		}
		for(int j = 0;j<kk;j++){
			int coff = 0;
			for(int l = 0;l<M;l++){
				coff += atom_new[l]*A_T[l*K+j];
			}
			for(int l = 0;l<M;l++){
				atom_new[l] -= coff*A_T[l*K+j];
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
		for(int i=0;i<M;i++){
			A_T[i*K+kk] = atom_new[i];
		}
		double* tmpA_T = new double[M*(kk+1)];
		for(int i =0 ;i<M;i++){
			for(int j = 0;j<kk+1;j++){
				tmpA_T[i*(kk+1)+j] = A_T[i*(K)+j];
			}
		}
		
		trmul(b,tmpA_T,1,M,kk+1,x_T);
		for(int i = 0;i<kk+1;i++){
			result[indx_set[i]] = x_T[i];
		}
		double *tmp = new double[M];
		//trmul(tmpA_T,x_T,M,kk+1,1,tmp);
		trmul(x_T,tmpA_T,1,kk+1,M,tmp);
		
		for(int i = 0;i<M;i++){
			r[i] = b[i] - tmp[i];
		}
		if (kk<K-1){
			trmul(r,A,1,M,N,Ar);
		}
	}
	double *pinv = new double[K*M];
   
	ginv(A_T_nonorth,pinv,M,K);
   
    trmul(b,pinv,1,K,M,x_T);
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
	double *omptmp = (double*)malloc(sizeof(double)*M);
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
        cout<<"2";
		OMP(y,Dtmp,M,K,int(*T),omptmp);
		for (int j = 0;j<K;j++){
			sc[i*K + j] = omptmp[j];
           // sc[i*K + j] =0;
		}
     
    }
	//delete []Dtmp;
	//delete []y;
    
}









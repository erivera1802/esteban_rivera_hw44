#include <stdlib.h>
#include <stdio.h>
#include <math.h>
float *ImportData(char *filename,int fil,int col);
float *MatrixMultiplication(float *A, float *B,int a_fil, int a_col,int b_fil,int b_col);
float *Transpose(float *A,int a_fil, int a_col);
float *Cholesky(float *A, int n);
float * SolverL(float*A, float*b,float*x,int n);
float * SolverU(float*C, float*d,float*y,int n);

int main(int argc, char **argv)
{
  
  
  //char filename[100]="new_data.dat";
  int i,j;
  float *YT;
  float *ThPh;
  float *G;
  float *Gt;
  float *A;
  float *Y;
  float *t;
  float *b;
  float *L;
  float *Lu;
  float *J;
  float *M;
  FILE *am;
  float *datas1;
  float *theta;
  char *datas="Theta_Phi.dat";
  char *punto1="punto1.dat";
  Y=malloc(sizeof(float)*38);
  t=malloc(sizeof(float)*38);
  theta=malloc(sizeof(float)*1000);
  Gt=malloc(sizeof(float)*38*2);
  A=malloc(sizeof(float)*3*3);
  b=malloc(sizeof(float)*3);
  L=malloc(sizeof(float)*3*3);
  Lu=malloc(sizeof(float)*3*3);
  J=malloc(sizeof(float)*3);
  M=malloc(sizeof(float)*3);
  am=malloc(sizeof(char)*1000);
  YT=ImportData(argv[1],39,2);
  
  datas1=ImportData(datas,1000,2);
  for(i=0;i<38;i++)
    {
      t[i]=YT[2*i];
      Y[i]=YT[2*i+1];
    }
  
  //printf("theta");
    for(i=0;i<100;i++)
    {
      theta[i]=datas1[2*i];
      //   printf("%f \n",theta[i]); 
    }
   
  for(i=0;i<38;i++)
    {
      for(j=0;j<3;j++)
	{
	  if(j==0)
	    {
	      G[i*3+j]=1;
	    }
	  else if(j==1)
	    {
	      G[i*3+j]=t[i];
	    }
	  else if(j==2)
	    {
	      G[i*3+j]=pow(t[i],2)/2;
	    }
	}
    }

  Gt=Transpose(G,38,3);
 
  A=MatrixMultiplication(Gt,G,3,38,38,3);

  b=MatrixMultiplication(Gt,Y,3,38,38,1);
 

  L=Cholesky(A,3);


  Lu=Transpose(L,3,3);

  J=SolverL(L,b,J,3);
  M=SolverU(Lu,J,M,3);
  printf("M\n");
    for (i=0;i<3;i++)
    { 
      printf("%f\n",M[i]);
    }    
  am=fopen(punto1,"a");
  if(!am)
    {
      printf("problem opening");
    }
    
  
  for(i=0;i<3;i++)
    {
      fprintf(am,"%f\t",M[i]);
    }
  fprintf(am,"\n");
  fclose(am);
  
    return 0;
}


float * SolverL(float*A, float*b,float*x,int n)
{
  int i;
  int j;
 for(i=0;i<n;i++)
    {
      x[i]=b[i];
      for(j=0;j<i;j++)
	{
	  x[i]-=A[i*n+j]*x[j];
	}
      x[i]=x[i]/A[i*n+i];
    }
 return x;
}

float * SolverU(float*A, float*b,float*x,int n)
{
  int i;
  int j;
 for(i=n-1;i>=0;i--)
    {
      x[i]=b[i];
      for(j=i+1;j<n;j++)
	{
	  x[i]-=A[i*n+j]*x[j];
	}
      x[i]=x[i]/A[i*n+i];
    }
 return x;
}



float *Cholesky(float *A,int n)
{
  int i;
  int j;
  int k;
  float suma1;
  float suma2;
   float *L;
   L=malloc(sizeof(float)*n*n);
  for (j=0;j<n;j++)
    {
      suma1=0;
      for(k=0;k<j;k++)
	{
	  suma1+=L[j*n+k]*L[j*n+k];
	}
      L[j*n+j]=sqrt(A[j*n+j]-suma1);
      for(i=j+1;i<n;i++)
	{
	  suma2=0;
	  for(k=0;k<j;k++)
	    {
	      suma2+=L[i*n+k]*L[j*n+k];
	    }
	  L[i*n+j]=(A[i*n+j]-suma2)/L[j*n+j]; 
	}
    }
  return L;
}



float *ImportData(char *filename ,int fil,int col)
{
  FILE *in;
  int i;
  int j=0;
  float var1,var2;
  int test;
  float *A;
  A=malloc(sizeof(float)*fil*col+2);

  in = fopen(filename, "r");
  if(!in)
    {
      printf("problems opening the file %s\n", filename);
      exit(1);
      
    }
  //  printf("Now I am reading\n");

  do
    {
      test = fscanf(in, "%f %f\n", &var1,&var2);
      if(test!=EOF)
	{
	  //	  printf("%f %f\n", var1,var2);
	  //printf("%d \n",j);
	  A[j]=var1;
	  A[j+1]=var2;
	  j=j+2;
	}
    }
  while(test!=EOF);
  fclose(in);
  // printf("%d \n",j);
  return A;
}

float *Transpose(float *A,int a_fil, int a_col)
{
  int i;
  int j;
  float *T;
  T=malloc(sizeof(float)*a_fil*a_col);
  for(i=0;i<a_fil;i++)
    {
      for(j=0;j<a_col;j++)
	{
	  T[j*a_fil+i]=A[i*a_col+j];
	}
    }
  return T;
}


float *MatrixMultiplication(float *A, float *B,int a_fil, int a_col,int b_fil,int b_col)
{
  int i;
  int j;
  int k;
  int m;
  float *C;
  float suma;
  C=malloc(sizeof(float)*a_fil*b_col);     
  if(a_col!=b_fil)
    {
      printf("Matrices cannot be multiplied");
      exit(1);
    }
  else
    {
      m=a_col;
    }
  for (i=0;i<a_fil;i++)
    {
      for(j=0;j<b_col;j++)
	{
	  suma=0;
	  for(k=0;k<m;k++)
	    {
	      suma+=A[i*a_col+k]*B[k*b_col+j];
	    }
	  C[i*b_col+j]=suma;
	}
    }
  return C;
}

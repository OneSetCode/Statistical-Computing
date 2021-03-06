#include "matrices.h"

//allocates the memory for a matrix with 
//n rows and p columns
double ** allocmatrix(int n,int p)
{
	int i;
	double** m;
	
	m = new double*[n];
	for(i=0;i<n;i++)
	{
		m[i] = new double[p];
		memset(m[i],0,p*sizeof(double));
	}
	return(m);
}

//frees the memory for a matrix with n rows
void freematrix(int n,double** m)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		delete[] m[i]; m[i] = NULL;
	}
	delete[] m; m = NULL;
	return;
}

//creates the copy of a matrix with n rows and p columns
void copymatrix(int n,int p,double** source,double** dest)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			dest[i][j] = source[i][j];
		}
	}
	return;
}

//reads from a file a matrix with n rows and p columns
void readmatrix(char* filename,int n,int p,double* m[])
{
	int i,j;
	double s;
	FILE* in = fopen(filename,"r");
	
	if(NULL==in)
	{
		printf("Cannot open input file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			fscanf(in,"%lf",&s);
			m[i][j] = s;
		}
	}
	fclose(in);
	return;
}

//prints the elements of a matrix in a file
void printmatrix(char* filename,int n,int p,double** m)
{
	int i,j;
	double s;
	FILE* out = fopen(filename,"w");
	
	if(NULL==out)
	{
		printf("Cannot open output file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		fprintf(out,"%.3lf",m[i][0]);
		for(j=1;j<p;j++)
		{
			fprintf(out,"\t%.3lf",m[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	return;
}

//creates the transpose of the matrix m
double** transposematrix(int n,int p,double** m)
{
	int i,j;
	
	double** tm = allocmatrix(p,n);
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<n;j++)
		{
			tm[i][j] = m[j][i];
		}
	}	
	
	return(tm);
}

//calculates the dot (element by element) product of two matrices m1 and m2
//with n rows and p columns; the result is saved in m
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			m[i][j] = m1[i][j]*m2[i][j];
		}
	}
	
	return;
}

//calculates the product of a nxp matrix m1 with a pxl matrix m2
//returns a nxl matrix m
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<n;i++)
	{
		for(k=0;k<l;k++)
		{
			s = 0;
			for(j=0;j<p;j++)
			{
				s += m1[i][j]*m2[j][k];
			}
			m[i][k] = s;
		}
	}
	return;
}

void set_mat_identity(int p, double *A)
{
 int i;

 for(i = 0; i < p * p; i++) A[i] = 0;
 for(i = 0; i < p; i++) A[i * p + i] = 1;
 return;
}

//computes the inverse of a symmetric positive definite matrix
void inverse(int p,double** m)
{
  int i,j,k;
  double* m_copy = (double*)malloc((p * p) * sizeof(double));
  double* m_inv = (double*)malloc((p * p) * sizeof(double));

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m_copy[k] = m[i][j];
        k++;
     }
  }

  set_mat_identity(p, m_inv);

  //-----  Use LAPACK  -------
  if(0!=(k=clapack_dposv(CblasRowMajor, CblasUpper, p, p, m_copy, p, m_inv, p)))
  {
    fprintf(stderr,"Something was wrong with clapack_dposv [%d]\n",k);
     exit(1);
  }
  //--------------------------

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m[i][j] = m_inv[k];
        k++;
     }
  }  

  free(m_copy);
  free(m_inv);

  return;
}


//computes the log of the determinant of a symmetric positive definite matrix
double logdet(int p,double** m)
{
	int i,j;
	char jobvl = 'N';
	char jobvr = 'N';
	int lda = p;
	double wr[2*p];
	double wi[2*p];
	double vl[p][p];
	int ldvl = p*p;
	double vr[p][p];
	int ldvr = p*p;
	double work[p*p];
	int lwork = p*p;
	double a[p][p];
	int info;
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dgeev_(&jobvl,&jobvr,&p,(double*)a,&lda,(double*)wr,(double*)wi,(double*)vl, 
		  &ldvl,(double*)vr,&ldvr,(double*)work,&lwork,&info);

	if(0!=info)
	{
		printf("Smth wrong in the call of 'dgeev' error is [info = %d]\n",info);
		exit(1);
	}	   
	
	double logdet = 0;
	for(i=0;i<p;i++) logdet+=log(wr[i]);	
	return(logdet);
}


//Implement the marglik function
double marglik(int n,int p,double** data,int lenA,int* A)
{
	double logml;

	double** D1 = allocmatrix(n, 1);
	double** DA = allocmatrix(n, lenA);
	double** MA = allocmatrix(lenA, lenA);

	for(int i = 0; i < n; i++)
	{
		D1[i][0] = data[i][0];
	}

	for(int i = 0; i < lenA; i++)
	{
		for(int j = 0; j < n; j++)
		{
			DA[j][i] = data[j][A[i] - 1];
		}
	}

	double** tDA = transposematrix(n, lenA, DA);

	matrixproduct(lenA, n, lenA, tDA, DA, MA);

	for(int i = 0; i < lenA; i++)
	{
		MA[i][i] += 1;
	}

	double** tD1 = transposematrix(n, 1, D1);
	double** tD1_D1 = allocmatrix(1, 1);
	matrixproduct(1, n, 1, tD1, D1, tD1_D1);

	double** tD1_DA = allocmatrix(1, lenA);
	matrixproduct(1, n, lenA, tD1, DA, tD1_DA);

	double** iMA = allocmatrix(lenA, lenA);
    copymatrix(lenA, lenA, MA, iMA);
	inverse(lenA, iMA);

	double** tDA_D1 = allocmatrix(lenA, 1);
	matrixproduct(lenA, n, 1, tDA, D1, tDA_D1);

	double** tD1_DA_iMA = allocmatrix(1, lenA);
	matrixproduct(1, lenA, lenA, tD1_DA, iMA, tD1_DA_iMA);

	double** tD1_DA_iMA_tDA_D1 = allocmatrix(1, 1);
	matrixproduct(1, lenA, 1, tD1_DA_iMA, tDA_D1, tD1_DA_iMA_tDA_D1);

	logml = (lgamma((n + lenA + 2.0) / 2.0) - lgamma((lenA + 2.0) / 2.0) - 
	    0.5 * logdet(lenA, MA) - ((n + lenA + 2.0) / 2.0) * 
		log(1.0 + **tD1_D1 - **tD1_DA_iMA_tDA_D1));
	
	freematrix(n,D1);
	freematrix(n,DA);
	freematrix(lenA, MA);
	freematrix(lenA, tDA);
	freematrix(1, tD1);
	freematrix(1, tD1_D1);
	freematrix(1, tD1_DA);
	freematrix(lenA, iMA);
	freematrix(lenA, tDA_D1);
	freematrix(1, tD1_DA_iMA);
	freematrix(1, tD1_DA_iMA_tDA_D1);

	return(logml);
}

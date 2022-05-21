#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>

gsl_matrix * subMatrix(gsl_matrix * M, int k) 
{
	int n = M -> size1;
	gsl_matrix * m = gsl_matrix_alloc(n - 1, n - 1);

	for(int i = 1; i < n; i++) 
	{
		for(int j = 0; j < n; j++) 
		{
			if(j < k)
			{
				gsl_matrix_set(m, i - 1, j, gsl_matrix_get(M, i, j));
			} 
			else 
			{
				gsl_matrix_set(m, i - 1, j - 1, gsl_matrix_get(M, i, j));
			}
		}
	}
	return(m);
}


double getDeterminant(gsl_matrix * M)
{
	int n = M -> size1;
	gsl_matrix * m = gsl_matrix_alloc(n-1, n-1);
	double a;
	double det = 0;

	if(n == 1) 
	{
		return(gsl_matrix_get(M, 0, 0));
	}

	if(n == 2)
	{
		return(gsl_matrix_get(M, 0, 0) * gsl_matrix_get(M, 1, 1) - gsl_matrix_get(M, 0, 1) * gsl_matrix_get(M, 1, 0));
	}

	for(int j = 0; j < n; j++) 
	{
		m = subMatrix(M, j);
		det += gsl_matrix_get(m, 0, j) * pow(-1, j + 2) * getDeterminant(m);
	}
}


int main()
{
	double det;

	gsl_matrix * M = gsl_matrix_alloc(10, 10);

	FILE * f = fopen("mybandedmatrix.txt", "r");
	gsl_matrix_fscanf(f, M);
	fclose(f);

	det = getDeterminant(M);

	gsl_matrix_free(M);
	printf("asfasfafa\n");

	printf("The determinant is %.4lf\n", det);
	return(det);
}

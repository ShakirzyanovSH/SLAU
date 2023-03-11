#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <stdio.h>

const
	int n=3;

void Gauss(double x[3], double e[3], double A[3][4])
{
	int i,j,k;
	double t;
	double A1[3][4];
	for(i=0;i<n;i++) for(j=0;j<=n;j++) A1[i][j] = A[i][j];
	
// проход в прямую сторону	
	for(i=0;i<n;i++)
	{
		for(j=i;j<n;j++)
		{
			if(j==i)
			{
				t = A1[j][j];
				for(k=0;k<=n;k++) A1[j][k] /= t;
			}
			else
			{
				t = A1[j][i];
				for(k=0;k<=n;k++) A1[j][k] -= t*A1[i][k];
			}
		}
	}
// проход в обратную сторону
	x[n-1] = A1[n-1][n];
	for(i=n-2;i>=0;i--)
	{
		t = 0;
		for(j=n-1;j>i;j--) t += x[j]*A1[i][j];
		x[i] = A1[i][n]-t;
	}
}

void Iteration(double x[3], double e[3], double A[3][4], double eIt)
{
	int i,j,k;
	double t;
	bool b=false;
	double A1[3][4];
	double xs[3];
	for(i=0;i<n;i++) for(j=0;j<=n;j++) A1[i][j] = A[i][j];

// начальное приближение
	for(i=0;i<n;i++) x[i] = 0;

	do
	{
		for(i=0;i<n;i++) xs[i] = x[i];
		for(i=0;i<n;i++)
		{
			t=0;
			for(j=0;j<n;j++)
			{
				if(j!=i) t += A1[i][j]*xs[j]/A1[i][i];
			}
			x[i] = A1[i][n]/A1[i][i] - t;	
	
		}
		printf(" x = (%-5.3lf %-5.3lf %-5.3lf) \n",x[0],x[1],x[2]);

		b=true;
		for(i=0;i<n;i++)
		{
			if(fabs(x[i]-xs[i])>eIt) b = false;
		}
	} while(!b);
}


void Zeidel(double x[3], double e[3], double A[3][4], double eIt)
{
	int i,j,k;
	double t;
	bool b=false;
	double A1[3][4];
	double xs[3], xprev[3];
	for(i=0;i<n;i++) for(j=0;j<=n;j++) A1[i][j] = A[i][j];

// начальное приближение
	for(i=0;i<n;i++) x[i] = 0;
	do
	{
		for(i=0;i<n;i++) 
		{
			xs[i] = x[i];
			xprev[i] = x[i];
		}
		for(i=0;i<n;i++)
		{
			t=0;
			for(j=0;j<n;j++)
			{
				if(j!=i) t += A1[i][j]*xs[j]/A1[i][i];
			}
			x[i] = A1[i][n]/A1[i][i] - t;	
			xs[i] = x[i];
		}
		printf(" x = (%-5.3lf %-5.3lf %-5.3lf) \n",x[0],x[1],x[2]);

		b=true;
		for(i=0;i<n;i++)
		{
			if(fabs(x[i]-xprev[i])>eIt) b = false;
		}
	} while(!b);
}


// проверка правильности решения
void RightSolve(double x[3], double e[3], double A[3][4])
{
	double t;	
	int i,j;
	for(i=0;i<n;i++)
	{
		t = 0;
		for(j=0;j<n;j++) t += x[j]*A[i][j];
// находим невязку		
		e[i] = A[i][n] - t;
	}
}

int main()
{

	int i,j;
	double A[n][n+1];
	double x[n];
	double e[n];
	double eIt;
	
	printf(" Solution of the system of linear equations. \n");
	printf(" 3.6x1 +2.4*x2+1.9x3 = 5 \n");
	printf(" 2.1x1 +7.5*x2+6.3x3 = 5.6 \n");
	printf(" 4.4x1 +5.1*x2+7.1x3 = 4.6 \n");

	FILE *fr = fopen("slau.txt","r");
	if(fr!=NULL)
	{
		for(i=0;i<n;i++)
		{
			for(j=0;j<=n;j++) fscanf(fr,"%lf",&A[i][j]); 
		}			
		
		printf("\n 1. Gauss method. \n");
		Gauss(x,e,A);
		printf(" System Solution: \n");
		printf(" x = \n");
		for(i=0;i<n;i++) printf(" %lf \n",x[i]);
		printf("\n Check the correctness of the found solution. \n");
		RightSolve(x,e,A);
		printf(" Discrepancy in solving the system: \n");
		printf(" e = \n");
		for(i=0;i<n;i++) printf(" %-7.2e \n",e[i]);


		printf("\n 2. Iterative Methods. \n");
		printf(" Specify Solution Accuracy for Iterative Methods e = ");
		scanf("%lf",&eIt);

		printf("\n 2.1. Simple iteration method. \n");
		Iteration(x,e,A,eIt);
		printf(" System Solution: \n");
		printf(" x = \n");
		for(i=0;i<n;i++) printf(" %lf \n",x[i]);
		printf("\n Check the correctness of the found solution. \n");
		RightSolve(x,e,A);
		printf(" Discrepancy in solving the system: \n");
		printf(" e = \n");
		for(i=0;i<n;i++) printf(" %-7.2e \n",e[i]);

		printf("\n 2.2. Seidel method. \n");
		Zeidel(x,e,A,eIt);
		printf(" System Solution: \n");
		printf(" x = \n");
		for(i=0;i<n;i++) printf(" %lf \n",x[i]);
		printf("\n Check the correctness of the found solution. \n");
		RightSolve(x,e,A);
		printf(" Discrepancy in solving the system: \n");
		printf(" e = \n");
		for(i=0;i<n;i++) printf(" %-7.2e \n",e[i]);
	}
	else
	{
		printf(" System data file not found \n");
	}

	getch();
	return 0;
}


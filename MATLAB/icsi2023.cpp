/*
  ICSI 2023 Special Session on Single Objective Bounded Optimization Problems

  Step 1: Run the following command in Matlab:
  mex icsi2023.cpp -DWINDOWS

  Step2: You can use the following statement in your program to use the test set:
  f = icsi2023(x, func_num);
  where x is a m*D matrix. D is the dimension of test functions. m is the
  number of tests at one time. func_num is from 1 to 10.
*/

#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

double *O, *M, *y, *z;

#define INF 1.0e99
#define EPS 1.0e-8
#define E 2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

// define head for icsi2022
// define for basic functions
void ellips_func(double *, double *, int, int, double *, int, double *);	  /* Ellipsoidal */
void bent_cigar_func(double *, double *, int, int, double *, int, double *);  /* Discus */
void discus_func(double *, double *, int, int, double *, int, double *);	  /* Bent_Cigar */
void rosenbrock_func(double *, double *, int, int, double *, int, double *);  /* Rosenbrock's */
void ackley_func(double *, double *, int, int, double *, int, double *);	  /* Ackley's */
void rastrigin_func(double *, double *, int, int, double *, int, double *);	  /* Rastrigin's  */
void griewank_func(double *, double *, int, int, double *, int, double *);	  /* Griewank's  */
void schwefel_func(double *, double *, int, int, double *, int, double *);	  /* Schwefel's */
void grie_rosen_func(double *, double *, int, int, double *, int, double *);  /* Griewank-Rosenbrock  */
void happycat_func(double *, double *, int, int, double *, int, double *);	  /* HappyCat */
void weierstrass_func(double *, double *, int, int, double *, int, double *); /* Weierstrass's  */
void katsuura_func(double *, double *, int, int, double *, int, double *);	  /* Katsuura */

void alple_func(double *, double *, int, int, double *, int, double *); /* alple 01 */
void dixon_func(double *, double *, int, int, double *, int, double *); /* Dixon-Prince's Function */

// not use

void ex3_func(double *, double *, int, int, double *, int, double *);
void logexp_func(double *, double *, int, int, double *, int, double *);
void invert_cos_wave_func(double *, double *, int, int, double *, int, double *);
void patho_func(double *, double *, int, int, double *, int, double *);
void salomon_func(double *, double *, int, int, double *, int, double *);
void sargan_func(double *, double *, int, int, double *, int, double *);
void wavy_func(double *, double *, int, int, double *, int, double *);
void hgbat_func(double *, double *, int , int ,double *, int ,double *); /* HGBat */


// for composition functions

void cf01(double *, double *, int, int, double *, int, double *); /* Composition Function 1 */
void cf02(double *, double *, int, int, double *, int, double *); /* Composition Function 2 */
void cf03(double *, double *, int, int, double *, int, double *); /* Composition Function 3 */
// void cf04 (double *, double *, int ,  int ,double *, int ,double *); /* Composition Function 4 */

// define for basic operator
void shiftfunc(double *, double *, int, double *);
void rotatefunc(double *, double *, int, double *);
void shift_rotate(double *, double *, int, double *, double *, double, int, int); /* shift and rotate */
// composition operator
void compose(double *, double *, double *, double *, int);

// end define

void icsi2023(double *, double *, int, int, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m, n, func_num;
	double *y, *x;
	// check the number of parameters
	if ((nrhs < 2) || (nlhs < 1))
	{
		mexErrMsgTxt("invalid number of parameters\n Example: y = icsi2022(x, func_num);\n");
	}
	// get data
	n = mxGetN(prhs[0]);
	m = mxGetM(prhs[0]);
	x = mxGetPr(prhs[0]);
	func_num = (int)*mxGetPr(prhs[1]);

	// check for dimension
	if (!(n == 10 || n == 20 || n == 50))
	{
		mexErrMsgTxt("Error: Test functions are only defined for D=10, 20, 50.");
	}
	
	// check for function index
	if (func_num > 10 || func_num < 1)
	{
		mexErrMsgTxt("Error: Test function is not defined.\n!");
	}
	// creat matrix for y
	plhs[0] = mxCreateDoubleMatrix(1, m, mxREAL);
	y = mxGetPr(plhs[0]);
	icsi2023(&x[0], &y[0], n, m, func_num);
}

void icsi2023(double *x, double *f, int nx, int mx, int func_num)
{
	int cf_num = 10, i, j;

	// init the data
	// ----------------------------------------
	FILE *filept;
	char FileName[256];
	free(M);
	free(O);
	free(y);
	free(z);
	y = (double *)malloc(sizeof(double) * nx);
	z = (double *)malloc(sizeof(double) * nx);

	/* Load Matrix M*/
	sprintf(FileName, "data/M_%d_D%d.txt", func_num, nx);
	filept = fopen(FileName, "r");
	if (filept == NULL)
	{
		printf("\nError: An error occurred while attempting open file %s!\n", FileName);
	}
	if (func_num < 8)
	{
		M = (double *)malloc(nx * nx * sizeof(double));
		if (M == NULL)
			printf("\nError: Insufficient Memory. Operation could not be completed!\n");
		for (i = 0; i < nx * nx; i++)
		{
			fscanf(filept, "%lf", &M[i]);
		}
	}
	else
	{
		M = (double *)malloc(cf_num * nx * nx * sizeof(double));
		if (M == NULL)
			printf("\nError: Insufficient Memory. Operation could not be completed!\n");
		for (i = 0; i < cf_num * nx * nx; i++)
		{
			fscanf(filept, "%lf", &M[i]);
		}
	}
	fclose(filept);

	/* Load shift_data */
	sprintf(FileName, "data/S_%d_D%d.txt", func_num, nx);
	filept = fopen(FileName, "r");
	if (filept == NULL)
	{
		printf("\nError: An error occurred while attempting open file %s!\n", FileName);
	}

	if (func_num < 8)
	{
		O = (double *)malloc(nx * sizeof(double));
		if (O == NULL)
			printf("\nError: Insufficient Memory. Operation could not be completed!\n");
		for (i = 0; i < nx; i++)
		{
			fscanf(filept, "%lf", &O[i]);
		}
	}
	else
	{
		O = (double *)malloc(nx * cf_num * sizeof(double));
		if (O == NULL)
			printf("\nError: Insufficient Memory. Operation could not be completed!\n");
		for (i = 0; i < cf_num - 1; i++)
		{
			for (j = 0; j < nx; j++)
			{
				fscanf(filept, "%lf", &O[i * nx + j]);
			}
			fscanf(filept, "%*[^\n]%*c");
		}
		for (j = 0; j < nx; j++)
		{
			fscanf(filept, "%lf", &O[nx * (cf_num - 1) + j]);
		}
	}
	fclose(filept);
	// ----------------------------------------
	// end init

	for (i = 0; i < mx; i++)
	{
		switch (func_num)
		{
		case 1:
			ellips_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=1000.0;
			break;
		case 2:
			bent_cigar_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=1000.0;
			break;
		case 3:
			rosenbrock_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=200.0;
			break;
		case 4:
			rastrigin_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=200.0;
			break;
		case 5:
			griewank_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=300;
			break;
		case 6:
			dixon_func(&x[i * nx], &f[i], nx, 1, O, 1,M);
			f[i]+=300;
			break;
		case 7:
			hgbat_func(&x[i * nx], &f[i], nx, 1,O, 0,M);
			f[i]+=500;
			break;
		case 8:
			cf01(&x[i * nx], &f[i], nx, 1,O, 1,M);
			break;
		case 9:
			cf02(&x[i * nx], &f[i], nx, 1,O, 1,M);
			break;
		case 10:
			cf03(&x[i * nx], &f[i], nx, 1,O, 1,M);
			break;

		default:
			f[i] = 0.0;
			break;
		}
	}
}

void bent_cigar_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Bent_Cigar */
{
	double index;
	int i;

	if(nx==50)
		index=6.0;
	else
		index=8.0;
	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag);
	f[0] = z[0] * z[0];
	for (i = 1; i < nx; i++)
	{
		f[0] += pow(10.0, index) * z[i] * z[i];
	}
}

void weierstrass_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Weierstrass's  */
{
	int i, j, k_max;
	double sum, sum2, a, b;
	a = 0.5;
	b = 3.0;
	k_max = 20;
	f[0] = 0.0;

	shift_rotate(x, z, nx, O, M, 0.5 / 100.0, s_flag, r_flag); 

	for (i = 0; i < nx; i++)
	{
		sum = 0.0;
		sum2 = 0.0;
		for (j = 0; j <= k_max; j++)
		{
			sum += pow(a, j) * cos(2.0 * PI * pow(b, j) * (z[i] + 0.5));
			sum2 += pow(a, j) * cos(2.0 * PI * pow(b, j) * 0.5);
		}
		f[0] += sum;
	}
	f[0] -= nx * sum2;
}

void ellips_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Ellipsoidal */
{
	int i;
	f[0] = 0.0;
	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 
	for (i = 0; i < nx; i++)
	{
		f[0] += pow(10.0, 8.0 * i / (nx - 1)) * z[i] * z[i];
	}
}


void rosenbrock_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Rosenbrock's */
{
	int i;
	double tmp1, tmp2;
	f[0] = 0.0;
	shift_rotate(x, z, nx, O, M, 2.048 / 100.0, s_flag, r_flag); 
	z[0] += 2.0;											  
	for (i = 0; i < nx - 1; i++)
	{
		z[i + 1] += 2.0; 
		tmp1 = z[i] * z[i] - 2 * z[i + 1];
		tmp2 = z[i] - 2.0;
		f[0] += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}
}

void schwefel_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Schwefel's  */
{
	int i;
	double tmp;
	f[0] = 0.0;

	shift_rotate(x, z, nx, O, M, 1000.0 / 100.0, s_flag, r_flag); 

	for (i = 0; i < nx; i++)
	{
		z[i] += 4.209687462275036e+002;
		if (z[i] > 500)
		{
			f[0] -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			f[0] += tmp * tmp / nx;
		}
		else if (z[i] < -500)
		{
			f[0] -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			f[0] += tmp * tmp / nx;
		}
		else
			f[0] -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	f[0] += 4.189828872724338e+002 * nx;
}

void rastrigin_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Rastrigin's  */
{
	int i;
	f[0] = 0.0;

	shift_rotate(x, z, nx, O, M, 5.12 / 100.0, s_flag, r_flag); 

	for (i = 0; i < nx; i++)
	{
		f[0] += (z[i] * z[i] * z[i] * z[i] - 20 * cos(2.0 * PI * z[i]) + 20.0);
	}
}

void shiftfunc(double *x, double *xshift, int nx, double *O)
{
	int i;
	for (i = 0; i < nx; i++)
	{
		xshift[i] = x[i] - O[i];
	}
}

void rotatefunc(double *x, double *y, int nx, double *M)
{
	int i, j;
	for (i = 0; i < nx; i++)
	{
		y[i] = 0;
		for (j = 0; j < nx; j++)
		{
			y[i] = y[i] + x[j] * M[i * nx + j];
		}
	}
}

void shift_rotate(double *x, double *sr_x, int nx, double *O, double *M, double sh_rate, int s_flag, int r_flag) 
{
	int i;
	if (s_flag == 1)
	{
		if (r_flag == 1)
		{
			shiftfunc(x, y, nx, O);
			for (i = 0; i < nx; i++) 
			{
				y[i] = y[i] * sh_rate;
			}
			rotatefunc(y, sr_x, nx, M);
		}
		else
		{
			shiftfunc(x, sr_x, nx, O);
			for (i = 0; i < nx; i++) 
			{
				sr_x[i] = sr_x[i] * sh_rate;
			}
		}
	}
	else
	{

		if (r_flag == 1)
		{
			for (i = 0; i < nx; i++) 
			{
				y[i] = x[i] * sh_rate;
			}
			rotatefunc(y, sr_x, nx, M);
		}
		else
			for (i = 0; i < nx; i++) 
			{
				sr_x[i] = x[i] * sh_rate;
			}
	}
}


void grie_rosen_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Griewank-Rosenbrock  */
{
	int i;
	double tmp, tmp1, tmp2;
	f[0] = 0.0;

	shift_rotate(x, z, nx, O, M, 5.0 / 100.0, s_flag, r_flag); 

	z[0] += 1.0; 
	for (i = 0; i < nx - 1; i++)
	{
		z[i + 1] += 1.0; 
		tmp1 = z[i] * z[i] - z[i + 1];
		tmp2 = z[i] - 1.0;
		tmp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
		f[0] += (tmp * tmp) / 4000.0 - cos(tmp) + 1.0;
	}
	tmp1 = z[nx - 1] * z[nx - 1] - z[0];
	tmp2 = z[nx - 1] - 1.0;
	tmp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	;
	f[0] += (tmp * tmp) / 4000.0 - cos(tmp) + 1.0;
}

void discus_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Discus */
{
	int i;
	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 
	f[0] = pow(10.0, 6.0) * z[0] * z[0];
	for (i = 1; i < nx; i++)
	{
		f[0] += z[i] * z[i];
	}
}


void ackley_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Ackley's  */
{
	int i;
	double sum1, sum2;
	sum1 = 0.0;
	sum2 = 0.0;

	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 

	for (i = 0; i < nx; i++)
	{
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * PI * z[i]);
	}
	sum1 = -0.2 * sqrt(sum1 / nx);
	sum2 /= nx;
	f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
}
void happycat_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* HappyCat,cec 2014*/

{
	int i;
	double beta, r2, sum_z;
	beta = 1.0 / 8.0;

	shift_rotate(x, z, nx, O, M, 5.0 / 100.0, s_flag, r_flag); 

	r2 = 0.0;
	sum_z = 0.0;
	for (i = 0; i < nx; i++)
	{
		z[i] = z[i] - 1.0; 
		r2 += z[i] * z[i];
		sum_z += z[i];
	}
	f[0] = pow(fabs(r2 - nx), 2 * beta) + (0.5 * r2 + sum_z) / nx + 0.5;
}
void hgbat_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* hgbat */

{
	int i;
	double beta, r2, sum_z;
	beta = 1.0 / 8.0;

	shift_rotate(x, z, nx, O, M, 5.0 / 100.0, s_flag, r_flag); 

	r2 = 0.0;
	sum_z = 0.0;
	for (i = 0; i < nx; i++)
	{
		z[i] = z[i] - 1.0; 
		r2 += z[i] * z[i];
		sum_z += z[i];
	}
	f[0] = pow(r2 * r2 - sum_z * sum_z, 2 * beta) + (0.5 * r2 + sum_z) / nx + 0.5;
}

void katsuura_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Katsuura  */
{
	int i, j;
	double tmp, tmp1, tmp2, tmp3;
	f[0] = 1.0;
	tmp3 = pow(1.0 * nx, 1.2);

	shift_rotate(x, z, nx, O, M, 5.0 / 100.0, s_flag, r_flag); 

	for (i = 0; i < nx; i++)
	{
		tmp = 0.0;
		for (j = 1; j <= 32; j++)
		{
			tmp1 = pow(2.0, j);
			tmp2 = tmp1 * z[i];
			tmp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
		}
		f[0] *= pow(1.0 + (i + 1) * tmp, 10.0 / tmp3);
	}
	tmp1 = 10.0 / nx / nx;
	f[0] = f[0] * tmp1 - tmp1;
}

void griewank_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Griewank's  */
{
	int i;
	double s, p;
	s = 0.0;
	p = 1.0;

	shift_rotate(x, z, nx, O, M, 600.0 / 100.0, s_flag, r_flag); 

	for (i = 0; i < nx; i++)
	{
		s += pow(10, 5 * i / nx) * z[i] * z[i];
		p *= cos(z[i] / sqrt(1.0 + i));
	}
	f[0] = 1.0 + s / 4000.0 - p;
}

void alple_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;

	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 

	f[0] = 0.0;
	for (i = 0; i < nx ; i++)
	{
		f[0] += fabs(z[i]*sin(z[i])+0.1*z[i]);
	}
}

void dixon_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;

	shift_rotate(x, z, nx, O, M, 10.0/100.0, s_flag, r_flag); 

	f[0] = (z[0]-1)*(z[0]-1);
	for (i = 1; i < nx ; i++)
	{
		f[0] += (i+1)*(2*z[i]*z[i]-z[i-1])*(2*z[i]*z[i]-z[i-1]);
	}
}

void ex3_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;
	double tmp;

	shift_rotate(x, z, nx, O, M, 10.0/100.0, s_flag, r_flag); 

	f[0]=0.0;
	for (i = 0; i < nx ; i++)
	{
		if(z[i]==0)
			tmp=0;
		else
			tmp=pow(z[i],6.0)*(sin(1/z[i])+2);
		
		f[0] += tmp;
	}
}

void logexp_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;
	double tmp;

	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 

	f[0]=0.0;
	tmp=0.0;
	for (i = 0; i < nx ; i++)
	{
		tmp+=z[i]*z[i];
		f[0] += tmp;
	}
	f[0]=log(2*exp(-1e-4*tmp));
}


void invert_cos_wave_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;
	double tmp;

	shift_rotate(x, z, nx, O, M, 5.0/100.0 , s_flag, r_flag); 

	f[0]=0.0;
	// for (i = 0; i < nx ; i++)
	// {
	// 	z[i]=z[i]/20.0;
	// }
	for (i = 0; i < nx-1 ; i++)
	{
		tmp+=z[i]*z[i]+z[i+1]*z[i+1]+0.5*z[i]*z[i+1];
		f[0] += -(exp(-tmp/8)*cos(4*sqrt(tmp)));
	}
	f[0]+=9;
}

void patho_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;
	double tmp1,tmp2;

	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 

	f[0]=0.0;

	for (i = 0; i < nx-1 ; i++)
	{
		tmp1=sqrt(100*z[i]*z[i]+z[i+1]*z[i+1]);
		tmp2=1+1e-3*(z[i]*z[i]-2*z[i]*z[i+1]+z[i+1]*z[i+1]);
		f[0] += 0.5+((pow(sin(tmp1),2)-0.5)/tmp2);
	}
	f[0]*=1e3;
}
void salomon_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i;
	double tmp;

	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 

	f[0]=0.0;
	tmp=0.0;
	for (i = 0; i < nx ; i++)
	{
		tmp=z[i]*z[i];
	}
	tmp=sqrt(tmp);
	f[0] = 1-cos(2*PI*tmp)+0.1*tmp;
	// f[0]*=1e3;
}
void sargan_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i,j;
	double *tmp;
	tmp = (double *)malloc(nx * sizeof(double));
	double sum=0.0;

	shift_rotate(x, z, nx, O, M, 1.0, s_flag, r_flag); 

	for(i=0;i<nx;i++)
	{
		sum+=z[i];
	}
	
	for(i=0;i<nx;i++)
	{
		tmp[i]=sum-z[i];
	}


	f[0]=0.0;
	for (i = 0; i < nx ; i++)
	{
		f[0]+=nx*(z[i]*z[i]+0.4*z[i]*tmp[i]);
	}
	// f[0]*=1e3;
}
void wavy_func(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M)
{
	int i,j;
	double *tmp;
	double k=10.0;
	tmp = (double *)malloc(nx * sizeof(double));
	double sum=0.0;

	shift_rotate(x, z, nx, O, M, PI/100.0, s_flag, r_flag); 
	f[0]=0.0;
	// for(i=0;i<nx;i++)
	// {
	// 	z[i]=z[i]*PI/100;
	// }
	for(i=0;i<nx;i++)
	{
		f[0]+=1-cos(k*z[i])*exp(-0.5*z[i]*z[i]);
	}
	f[0]*=100.0;
	
}


// for caculate the composition function

void compose(double *x, double *f, double *bias, double *fit, int f_num)
{
	int i, j;
	double *w;
	double w_max = 0, w_sum = 0;
	double sigma = 2;
	w = (double *)malloc(f_num * sizeof(double));
	f[0] = 1;
	for (i = 0; i < f_num; i++)
	{
		fit[i] += bias[i];
		f[0] = f[0] * fit[i];
	}

	f[0] = 1e-3 * (exp(sigma * log(f[0] + 1)) - 1);
	free(w);
}

void cf01(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Composition Function 1 */
{
	int i, f_num = 4;
	double fit[4];
	double bias[4] = {0, 1e3, 1e2, 1e2};

	i = 0;
	ackley_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e4;
	i = 1;
	ellips_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e+10;
	i = 2;
	griewank_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e4;
	i = 3;
	rastrigin_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 3e1;
	compose(x, f, bias, fit, f_num);
}

void cf02(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Composition Function 2 */
{
	int i, f_num = 4;
	double fit[4];
	double bias[4] = {0, 1e1, 1e2, 1e3};
	i = 0;
	alple_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e3;
	i = 1;
	katsuura_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = 1*fit[i] / 5e1;
	i = 2;
	rosenbrock_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e1;
	i = 3;
	rastrigin_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e3;
	compose(x, f, bias, fit, f_num);
}

void cf03(double *x, double *f, int nx, int s_flag, double *O, int r_flag, double *M) /* Composition Function 3 */
{
	int i, f_num = 5;
	double fit[5];
	double bias[5] = {0, 1e2, 1e2, 1e2, 1e2};
	i = 0;
	happycat_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e+5;
	i = 1;
	grie_rosen_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e+5;
	i = 2;
	schwefel_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e+5;
	i = 3;
	weierstrass_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e+8;
	i = 4;
	ellips_func(x, &fit[i], nx, s_flag ,&O[i * nx],r_flag ,&M[i * nx * nx]);
	fit[i] = fit[i] / 1e+10;
	compose(x, f, bias, fit, f_num);
}


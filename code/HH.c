/* Simulation of Arnold tongues in RSA pacing 

TODO:

Measure time diff from swapping current array to current function
Check similarity between results with python script
Save current array out
*/

// #include "stdafx.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#define PI 3.14159265358979323846
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define MAXSTP 1e7
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double time_spent = 0.0;
double I_in(double, int );
void update_TDE_var();
double randn(double, double);
void derivs(double, double *, double *);
void rkck(double *, double *, double, double, double *, double *);
void rkqs(double *, double *, double *, double, double, double *, double *, double *);
void odeint(double *, double, double, double, double, double, short);
double hs(double);
double im(double), ih(double), in(double);
double tde_var_means[3];
double il, igm, igh, ign, alpha, itaum, itauh, itaun, iTm, iTh, iTn, vtm, vth, vtn, El, Iapp, Idark, cycle_start_t, sigma;
double betal, betam, betah, betan, beta, betataum, betatauh, betataun;
double ys[4];
double Tmod,dutycycle,rel_rsa,Iinj;
int neqs, tde_var_index;
int kount, kmax;
double maxarg1, maxarg2, minarg1, minarg2, dxsav, tstep;
FILE *pf, *pf2;

int main(int argc, char *argv[])
{
	// change folder name
	// figure out f vs I curve
	// Investigate params
	// decide of file dir struct
	// figure out if relaxation makes better results 
	// figure out hw to save TDE data wo seg error

	
	// main params: output_dir, sigma, tdeidx, i_base, rel_rsa, dutycycle, t_resp
	double hstep, hmin, eps, tstart, tend,TModStep,T0;
	int sflag, NTMod, i, j;
	clock_t begin = clock();

	// EDIT THIS
	
	/* Generating the current profile */
	
	neqs = 4;			/* Number of equations in the model */
	il = 0.23105;		/* Leak current (nA) (all capacitances are 1pF) */
	igm = 162.75;		/* Na activation current (nA) */
	igh = 4.638;		/* Na inactivation current (nA) */
	ign = 164.18;		/* K activation current (nA) */
	betal = 1.;			/* Slope of leak current (V^-1) */
	betam = 8.405;		/* Slope of Na activation */
	betah = 3.581;		/* Slope of Na inactivation */
	betan = 8.372;		/* slope of K activation */
	beta = 14;			/* VLSI transistor slope (V^-1) */
	betataum = 0.0;		/* Half width of Na bell-shape activation kinetics */
	betatauh = 0.0;		/* Half width of Na bell-shape inactivation kinetics */
	betataun = 0.0;		/* Half width of K bell-shape actiation kinetics */
	alpha = 39.54;		/* Coefficient of amplification of injected current */
	itaum = 0.6854;		/* Na activation time constant (nA) */
	itauh = 0.1482;		/* Na inactivation time constant (nA) */
	itaun = 0.6747;		/* K activation time constant (nA) */
	iTm = 0.;			/* Na activation bell shape amlitude (nA) */
	iTh = 0.;			/* Na inactivation bell shape amplitude (nA) */
	iTn = 0.;			/* K activation bell shape amplitude (nA) */
	vtm = 0.908;		/* Na activation threshold (V) */
	vth = 1.143;		/* Na inactivation threshold (V) */
	vtn = 0.911;		/* K activation threshold (V) */
	El = 0.6194;		/* Leak reversal potential */
	Idark = 0;			/* Dark current (nA) */
		
	Iinj = strtof(argv[4],NULL); /* Current injection during expiratory phase (nA) 0.0024 default */
	Tmod = strtof(argv[7],NULL);		/* Period of respiratory modulation (ms) 10 default */
	dutycycle = strtof(argv[6],NULL);	/* Duty cycle of respiration */
	rel_rsa = strtof(argv[5],NULL);			/* Relative amplitude of RSA (Imax - Iinj) / Iinj, where Imax is inhalation current */	
	tde_var_index = strtol(argv[3],NULL,10); /* specifies variable  to be set to vary wrt to time: 0 : Tmod, 1: dutycycle, 2: rsa */
	// step_var_index = strtof(argv[4],NULL);
	sigma = strtof(argv[2], NULL);         /* standard dev of variation for random variable tde_var */
	tde_var_means[0] = Tmod; /* defines mean values for variables */
	tde_var_means[1] = dutycycle;
	tde_var_means[2] = rel_rsa;
	
	sflag = 1;		/* save data */
	tstart = 0;
	tend = 2000;
	hstep = 0.01;			/* Trial step for integration - ms */
	hmin = 0.0;				/* Smallest value of integration step */
	dxsav = 0.02;			/* Maximum time interval between saves - ms */
	eps = 1e-12;			/* Accuracy on solution */
	kmax = 200000;			/* Maximum number of steps that can be stored */
	NTMod = 1000;
	T0 = 40;
	TModStep = (T0 - T0/16)/ (NTMod - 1);
	
	
	for (i=0;i<NTMod;i++)
		{	
			printf("start of loop\n");

			clock_t begin = clock();


			char outfilename[256];
			char outfilename2[256];

			char outfile[1024];
			char outfile2[1024];

			Tmod = T0/16 + i*TModStep;
			cycle_start_t = 0.; /* start time for injected current signal */
			// // generate file paths
			snprintf(outfilename, sizeof outfilename, "/sigma_%lf_tdeIdx_%i_iBase_%lf_relRsa_%lf_dc_%lf_tResp_%lf_tend_%lf.csv", sigma, tde_var_index, Iinj, rel_rsa, dutycycle, Tmod, tend);				 
			strcpy(outfile,argv[1]);
			strcat(outfile,outfilename);
			
			/*seg fault somewhere here*/
			// snprintf(outfilename2, sizeof outfilename2,"/TDE_var_tResp_%lf.csv",Tmod);
			// strcpy(outfile2,argv[1]);
			// strcat(outfile2,outfilename2);
		
			// /* Opening the output data file */
			if (sflag == 1)
			{
				pf = fopen(outfile, "w+");
				if (sigma !=0) pf2 = fopen(outfile2,"w+");
			}
			
			/* Generates the starting values ys by relaxation with constant current injection Iapp */
			for (j = 0; j < neqs; j++)  ys[j] = 0.466; //reset ys between run
			// sflag = 0;		/* do not save data */
			// odeint(ys, tstart, Tmod*100, eps, hstep, hmin, sflag); 
			
			/* Integrate equations */
			odeint(ys, tstart, tend, eps, hstep, hmin, sflag);
			clock_t end = clock();
			time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
			printf("Finished odeint loop %d, saveflag = %hu: saving to %s\nThe elapsed time is %f seconds\n\n",i, sflag, outfile, time_spent);

			fclose(pf);
			// fclose(pf2);

 
    		// printf("The elapsed time is %f seconds\n \n", time_spent);


		}
	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    // printf("The elapsed time is %f seconds\n", time_spent);
	
	return(0);
}

double I_in(double x,int saveflag)
	{
		if ( (x - cycle_start_t > Tmod) && sigma != 0 ) /* if new cycle and sigma set to non zero value, modify value for tde_var */
			{
				cycle_start_t = x;
				update_TDE_var(saveflag);
			}
		if (fmod((x-cycle_start_t)/Tmod,1.) <= dutycycle)
			{	
				return Iinj*(1+rel_rsa); /* within inhalation part of respiratory signal */
			}
		
		else 
			{ 
				return Iinj; /* within exhalatory part of respiratory signal */
			}
	}

void update_TDE_var(int saveflag)
	{
		double updated_val;
		updated_val = randn(tde_var_means[tde_var_index],sigma); /* generate a new value for tde_var based on a normal distribution */
		if (updated_val < 0) updated_val = 0;
		if (saveflag && sigma !=0) fprintf(pf2,"%lf\n",updated_val);
		if (tde_var_index == 0) Tmod = updated_val;
		else if (tde_var_index == 1) dutycycle = updated_val;
		else rel_rsa = updated_val;
	} 

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

double hs(double xx)
{
	double result;
	result = 0.5*(1 + tanh(100*xx));
	return(result);
}

double im(double vvm)
{
	double iim;

	iim = 0.5*igm*(1 + tanh(betam*(vvm - vtm)));

	return(iim);
}

double ih(double vvh)
{
	double iih;

	iih = 0.5*igh*(1 + tanh(betah*(vvh - vth)));

	return(iih);
}

double in(double vvn)
{
	double iin;

	iin = 0.5*ign*(1 + tanh(betan*(vvn - vtn)));

	return(iin);
}

void derivs(double z, double *yy, double *dyydx)
{
	double imyy1, ihyy2;

	imyy1 = im(yy[1]);
	ihyy2 = ih(yy[2]);
	dyydx[0] = (imyy1 - ihyy2)*hs(imyy1 - ihyy2) - in(yy[3]) + il * tanh(betal*(El - yy[0])) + alpha * Iapp + Idark;
	dyydx[1] = itaum * tanh(beta*(yy[0] - yy[1])) / (1 + (iTm / (4 * itaum))*(1-tanh(betataum*(yy[0] - vtm))*tanh(betataum*(yy[0] - vtm))));
	dyydx[2] = itauh * tanh(beta*(yy[0] - yy[2])) / (1 + (iTh / (4 * itauh))*(1-tanh(betatauh*(yy[0] - vth))*tanh(betatauh*(yy[0] - vth))));
	dyydx[3] = itaun * tanh(beta*(yy[0] - yy[3])) / (1 + (iTn / (4 * itaun))*(1-tanh(betataun*(yy[0] - vtn))*tanh(betataun*(yy[0] - vtn))));
}


void odeint(double *ystart, double x1, double x2, double eps, double h1, double hmin, short saveflag)
{
	double xsav, x, hnext, hdid, h, aa, bb;
	double *yscal, *y, *dydx;
	long nstp;
	short i;

	yscal = (double *)malloc(neqs * sizeof(double));
	y = (double *)malloc(neqs * sizeof(double));
	dydx = (double *)malloc(neqs * sizeof(double));

	x = x1;
	h = SIGN(h1, x2 - x1); 
	kount = 0;
	for (i = 0; i<neqs; i++)   y[i] = ystart[i];
	if (kmax>0)  xsav = x - dxsav * 2.0;			/* step back to store the first step */
	for (nstp = 0; nstp<MAXSTP; nstp++)
	{
		if (saveflag == 1)
		{												/* Current at time x */
			Iapp = I_in(x,saveflag);
		}
		else
		{
			Iapp = I_in(0,saveflag);						/* Initial current used in the relaxation step */
		}

		derivs(x, y, dydx);
		yscal[0] = 1.;								/* Amplitude of voltage oscillations */
		for (i = 1; i<neqs; i++)  yscal[i] = 1.;	/* Amplitude of gate parameter oscillations */
		if (kmax>0 && kount<kmax - 1 && fabs(x - xsav)>fabs(dxsav) && saveflag == 1)
		{
			fprintf(pf, "%lf,%lf,%lf,%lf,%lf,%lf\n", x, Iapp, y[0], y[1], y[2], y[3]);
			xsav = x;
			kount++;
		}

		if ((x + h - x2)*(x + h - x1)>0.0)   h = x2 - x;		/* if step size overshoots, decrease */
		rkqs(y, dydx, &x, h, eps, yscal, &hdid, &hnext);
		if ((x - x2)*(x2 - x1) >= 0.0)				/* Finished ? */
		{
			for (i = 0; i<neqs; i++)  ystart[i] = y[i];
			if (kmax && saveflag == 1)
			{
				fprintf(pf, "%lf,%lf,%lf,%lf,%lf,%lf\n", x, Iapp, y[0], y[1], y[2], y[3]);
				kount++;
			}

			free(dydx);
			free(y);
			free(yscal);
			return;
		}
		if (fabs(hnext) <= hmin)
		{
			printf("%s\n", "Step size too small in odeint");
			exit(-1);
		}
		h = hnext;
	}
	printf("%s\n", "Too many steps in routine odeint");
	exit(-1);
}



void rkqs(double *yy, double *dyydx, double *x, double htry, double eps, double *yyscal, double *hdid, double *hnext)
{
	short i;
	double errmax, h, htemp, xnew, *yyerr, *yytemp;

	yyerr = (double *)malloc(neqs * sizeof(double));
	yytemp = (double *)malloc(neqs * sizeof(double));
	h = htry;

	for (;;)
	{
		rkck(yy, dyydx, *x, h, yytemp, yyerr);
		errmax = 0.0;
		for (i = 0; i<neqs; i++) errmax = FMAX(errmax, fabs(yyerr[i] / yyscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp = SAFETY * h*pow(errmax, PSHRNK);
		h = (h >= 0.0 ? FMAX(htemp, 0.1*h) : FMIN(htemp, 0.1*h));
		xnew = (*x) + h;
		if (xnew == *x)
		{
			printf("%s\n", "stepsize underflow in rkqs");
			exit(-1);
		}
	}
	if (errmax > ERRCON) *hnext = SAFETY * h*pow(errmax, PGROW);
	else *hnext = 5.0*h;
	*x += (*hdid = h);
	for (i = 0; i<neqs; i++) yy[i] = yytemp[i];

	free(yytemp);
	free(yyerr);
}


void rkck(double *yy, double *dyydx, double x, double hh, double *yyout, double *yyerr)
{
	short i;
	static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0, dc5 = -277.0 / 14336.0;
	double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0, dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;
	double *ak2, *ak3, *ak4, *ak5, *ak6, *yytemp;

	ak2 = (double *)malloc(neqs * sizeof(double));
	ak3 = (double *)malloc(neqs * sizeof(double));
	ak4 = (double *)malloc(neqs * sizeof(double));
	ak5 = (double *)malloc(neqs * sizeof(double));
	ak6 = (double *)malloc(neqs * sizeof(double));
	yytemp = (double *)malloc(neqs * sizeof(double));

	for (i = 0; i<neqs; i++)   yytemp[i] = yy[i] + b21 * hh*dyydx[i];
	derivs(x + a2 * hh, yytemp, ak2);
	for (i = 0; i<neqs; i++)   yytemp[i] = yy[i] + hh * (b31*dyydx[i] + b32 * ak2[i]);
	derivs(x + a3 * hh, yytemp, ak3);
	for (i = 0; i<neqs; i++)   yytemp[i] = yy[i] + hh * (b41*dyydx[i] + b42 * ak2[i] + b43 * ak3[i]);
	derivs(x + a4 * hh, yytemp, ak4);
	for (i = 0; i<neqs; i++)   yytemp[i] = yy[i] + hh * (b51*dyydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
	derivs(x + a5 * hh, yytemp, ak5);
	for (i = 0; i<neqs; i++)   yytemp[i] = yy[i] + hh * (b61*dyydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
	derivs(x + a6 * hh, yytemp, ak6);
	for (i = 0; i<neqs; i++)   yyout[i] = yy[i] + hh * (c1*dyydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
	for (i = 0; i<neqs; i++)   yyerr[i] = hh * (dc1*dyydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);

	free(yytemp);
	free(ak6);
	free(ak5);
	free(ak4);
	free(ak3);
	free(ak2);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void evolve_inhom_constrained(int n_x, int n_y, double delta_x,
							  double delta_y, double kappa, double A, double delta_t, int time_steps,
							  double t0, fftw_complex *comp, double Ceff[3][3][3][3],
							  double DeltaC[3][3][3][3], double **sigma_T, double **epsilon_T,
							  double ave_comp, double **sig_app, int STEPS, int MAXITR, double MAXERR,
							  int n_alpha, int n_beta, double *alpha, double *beta, double *alpha_prime,
							  double *beta_prime, int plan_indicator, int SNA, int NOISEADDER,
							  double noise_str)
{

	FILE *fpd;
	FILE *fptmp;
	FILE *fpstr;

	int INDEX = 0;
	int i1, i2, i3, i4, i5, i6;
	int J;

	int half_nx, half_ny;
	double *kx, *ky;
	double delta_kx, delta_ky;
	double k2, k4;
	double inv_denom;
	double gamma_I = 25.0;
	double gamma_A = -20.0;

	int nxny;
	double inv_nxny;

	char NAME[50];

	double epsstr11, epsstr22, epsstr12;

	double *c;
	double **Del_sigma_T;
	double **E;
	double **strain;
	double *Omega11;
	double *Omega12;
	double *Omega21;
	double *Omega22;
	double S[3][3][3][3];

	double temp1, temp2, temp3;

	double a, b;
	double ap, bp;
	double realc;
	size_t tmp;

	unsigned FLAG;

	fftw_complex *g;
	fftw_complex *u1_old, *u2_old;
	fftw_complex *u1_new, *u2_new;
	fftw_complex *eps_star11, *eps_star12, *eps_star22;
	fftw_complex *eta;
	fftw_complex *mu_el;

	fftw_plan planF, planB;

	nxny = n_x * n_y;
	inv_nxny = 1.0 / nxny;

	g = fftw_malloc(nxny * sizeof(fftw_complex));
	u1_old = fftw_malloc(nxny * sizeof(fftw_complex));
	u2_old = fftw_malloc(nxny * sizeof(fftw_complex));
	u1_new = fftw_malloc(nxny * sizeof(fftw_complex));
	u2_new = fftw_malloc(nxny * sizeof(fftw_complex));
	eps_star11 = fftw_malloc(nxny * sizeof(fftw_complex));
	eps_star12 = fftw_malloc(nxny * sizeof(fftw_complex));
	eps_star22 = fftw_malloc(nxny * sizeof(fftw_complex));
	eta = fftw_malloc(2 * 2 * sizeof(fftw_complex));
	mu_el = fftw_malloc(nxny * sizeof(fftw_complex));

	if (plan_indicator == 0)
	{
		FLAG = FFTW_ESTIMATE;
	}
	else if (plan_indicator == 1)
	{
		FLAG = FFTW_MEASURE;
	}
	else if (plan_indicator == 2)
	{
		FLAG = FFTW_PATIENT;
	}
	else if (plan_indicator == 3)
	{
		FLAG = FFTW_EXHAUSTIVE;
	}

	// fftw_plan_with_nthreads(8);

	planF =
		fftw_plan_dft_2d(n_x, n_y, mu_el, mu_el, FFTW_FORWARD, FLAG);
	planB =
		fftw_plan_dft_2d(n_x, n_y, eps_star11, eps_star11, FFTW_BACKWARD, FLAG);

	/* Declare and initialise the acoustic tensor */

	half_nx = (int)n_x / 2;
	half_ny = (int)n_y / 2;

	kx = (double *)malloc((size_t)n_x * sizeof(double));
	ky = (double *)malloc((size_t)n_y * sizeof(double));

	delta_kx = (2.0 * M_PI) / (n_x * delta_x);
	delta_ky = (2.0 * M_PI) / (n_y * delta_y);

	Omega11 = dvector(0, nxny);
	Omega12 = dvector(0, nxny);
	Omega21 = dvector(0, nxny);
	Omega22 = dvector(0, nxny);

	for (i1 = 0; i1 < n_x; ++i1)
	{
		if (i1 < half_nx)
			kx[i1] = i1 * delta_kx;
		else
			kx[i1] = (i1 - n_x) * delta_kx;
	}
	for (i2 = 0; i2 < n_y; ++i2)
	{
		if (i2 < half_ny)
			ky[i2] = i2 * delta_ky;
		else
			ky[i2] = (i2 - n_y) * delta_ky;
	}

	for (i1 = 0; i1 < n_x; ++i1)
	{
		for (i2 = 0; i2 < n_y; ++i2)
		{
			J = i2 + n_y * i1;
			Omega11[J] = 0.0;
			Omega12[J] = 0.0;
			Omega21[J] = 0.0;
			Omega22[J] = 0.0;
		}
	}

	/* Write the initial configuration */

	c = (double *)malloc((size_t)nxny * sizeof(double));
	for (i1 = 0; i1 < n_x; ++i1)
	{
		for (i2 = 0; i2 < n_y; ++i2)
		{
			J = i2 + n_y * i1;
			c[J] = creal(comp[J]);
		}
	}

	sprintf(NAME, "../output/data/time%d.dat",
			(int)(INDEX + t0));
	fpd = fopen(NAME, "w");
	tmp = fwrite(&c[0], sizeof(double), (size_t)nxny, fpd);
	fclose(fpd);

	/** Calculate the Del_sigma_T tensor **/

	Del_sigma_T = dmatrix(1, 2, 1, 2);
	calculate_Del_sigma_T(DeltaC, epsilon_T, Del_sigma_T);
	E = dmatrix(1, 2, 1, 2);
	strain = dmatrix(1, 2, 1, 2);

	/* Calculate the homogeneous strain (if there is applied stress). In
	 * the absence of applied strain, we assume the homogeneous strain to
	 * be zero */

	calculate_S_exp(n_x, n_y, delta_x, delta_y, comp, ave_comp, Ceff, DeltaC, S, n_alpha,
					alpha);

	for (i1 = 1; i1 < 3; ++i1)
	{
		for (i2 = 1; i2 < 3; ++i2)
		{
			E[i1][i2] = 0.0;
		}
	}

	for (i1 = 1; i1 < 3; ++i1)
	{
		for (i2 = 1; i2 < 3; ++i2)
		{
			for (i3 = 1; i3 < 3; ++i3)
			{
				for (i4 = 1; i4 < 3; ++i4)
				{
					E[i1][i2] = E[i1][i2] + S[i1][i2][i3][i4] * sig_app[i3][i4];
				}
			}
		}
	}

	/*
	E[1][1] = E[2][2] = 0.01;
	E[1][2] = E[2][1] = 0.0;
	*/

	printf("%le\n", E[1][1]);
	printf("%le\n", E[2][2]);
	printf("%le\n", E[1][2]);
	printf("%le\n", E[2][1]);

	/* Evolve */

	fptmp = fopen("../output/noise_info", "a");

	for (INDEX = 0; INDEX < time_steps + 1; ++INDEX)
	{

		/* Calculate g and its Fourier transform */

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				J = i2 + n_y * i1;
				realc = creal(comp[J]);
				g[J] = 2.0 * A * realc * (1.0 - realc) * (1.0 - 2.0 * realc) + 0.0 * _Complex_I;
			}
		}

		fftw_execute_dft(planF, g, g);

		if (INDEX == 0)
		{

			/* The acoustic tensor Omega and the (homogeneous) base solution */

			calculate_Omega(n_x, n_y, half_nx, half_ny, kx, ky, Ceff,
							Omega11, Omega12, Omega21, Omega22);

			calculate_uzero(n_x, n_y, half_nx, half_ny, delta_kx, delta_ky, ave_comp,
							comp, Ceff, sigma_T, u1_old, u2_old, n_beta, beta, Omega11, Omega12, Omega21,
							Omega22, planF);
		}

		/* Refine the displacements */

		refine_u_chen(n_x, n_y, half_nx, half_ny, delta_x, delta_y, epsilon_T,
					  delta_kx, delta_ky, comp, ave_comp, Ceff, DeltaC, sigma_T,
					  u1_old, u2_old, u1_new, u2_new, Del_sigma_T, sig_app, E, MAXITR, MAXERR,
					  n_alpha, n_beta, alpha, beta, Omega11, Omega12, Omega21, Omega22, planF, planB);

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				J = i2 + n_y * i1;
				eps_star11[J] = _Complex_I * kx[i1] * u1_new[J];
				eps_star22[J] = _Complex_I * ky[i2] * u2_new[J];
				eps_star12[J] = 0.5 * _Complex_I * (kx[i1] * u2_new[J] + ky[i2] * u1_new[J]);
			}
		}

		/* Get the strain back to the real space */

		fftw_execute(planB);
		fftw_execute_dft(planB, eps_star12, eps_star12);
		fftw_execute_dft(planB, eps_star22, eps_star22);

		/* Calculate mu_el and take it to the Fourier space */
		if (INDEX != 0 && INDEX % STEPS == 0)
		{
			sprintf(NAME, "../output/data/stress-at-time%d.dat",
					(int)(INDEX + t0));
			fpstr = fopen(NAME, "w");
		}

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{

				J = i2 + n_y * i1;

				epsstr11 = creal(eps_star11[J] * inv_nxny);
				epsstr22 = creal(eps_star22[J] * inv_nxny);
				epsstr12 = creal(eps_star12[J] * inv_nxny);

				realc = creal(comp[J]);
				mu_el[J] = 0.0;

				a = alpha[(int)(n_alpha * realc)];
				b = beta[(int)(n_beta * realc)];

				ap = alpha_prime[(int)(n_alpha * realc)];
				bp = beta_prime[(int)(n_beta * realc)];

				strain[1][1] =
					E[1][1] + epsstr11 - epsilon_T[1][1] * b;
				strain[1][2] =
					E[1][2] + epsstr12 - epsilon_T[1][2] * b;
				strain[2][1] =
					E[2][1] + epsstr12 - epsilon_T[2][1] * b;
				strain[2][2] =
					E[2][2] + epsstr22 - epsilon_T[2][2] * b;

				temp1 = 0.;
				temp2 = 0.;
				temp3 = 0.;
				for (i3 = 1; i3 < 3; ++i3)
				{
					for (i4 = 1; i4 < 3; ++i4)
					{
						for (i5 = 1; i5 < 3; ++i5)
						{
							for (i6 = 1; i6 < 3; ++i6)
							{
								if (i3 == 1 && i4 == 1)
									temp1 = temp1 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
								if (i3 == 2 && i4 == 2)
									temp2 = temp2 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
								if ((i3 == 1 && i4 == 2) || (i3 == 2 && i4 == 1))
									temp3 = temp3 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
								mu_el[J] = mu_el[J] + 0.5 * ap * DeltaC[i3][i4][i5][i6] * strain[i5][i6] * strain[i3][i4] - bp * epsilon_T[i3][i4] * strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
							}
						}
					}
				}
				if (INDEX != 0 && INDEX % STEPS == 0)
				{
					fprintf(fpstr, "%d %d %le %le %le %le\n", i1, i2, realc, temp1, temp2, temp3);
				}
			}
		}
		if (INDEX != 0 && INDEX % STEPS == 0)
		{
			fclose(fpstr);
		}
		fftw_execute(planF);
		fftw_execute_dft(planF, comp, comp);

		/* Evolve the composition profile */

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{

				J = i2 + n_y * i1;

				k2 = kx[i1] * kx[i1] + ky[i2] * ky[i2];
				k4 = k2 * k2;

				inv_denom = 1.0 + 2. * k2 * delta_t * (kappa * k2 + gamma_I * k4 + gamma_A * (kx[i1] * kx[i1] * kx[i1] * kx[i1] + ky[i2] * ky[i2] * ky[i2] * ky[i2]));
				inv_denom = 1.0 / inv_denom;

				comp[J] = inv_denom * (comp[J] - k2 * delta_t * (g[J] + mu_el[J]));
			}
		}

		/* Get the composition back to real space */

		fftw_execute_dft(planB, comp, comp);
		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				J = i2 + n_y * i1;
				comp[J] = comp[J] * inv_nxny;
			}
		}

		/* Once in a while check that the composition is within bounds,
		 * and write the output **/

		if (INDEX != 0 && INDEX % STEPS == 0)
		{
			for (i1 = 0; i1 < n_x; ++i1)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					J = i2 + n_y * i1;
					c[J] = creal(comp[J]);
					if (c[J] <= -1.0 || c[J] >= 2.0)
					{
						printf("The composition goes out of bounds\n");
						printf("Exiting\n");
						exit(0);
					}
				}
			}

			sprintf(NAME, "../output/data/time%d.dat",
					(int)(INDEX + t0));
			fpd = fopen(NAME, "w");
			tmp = fwrite(&c[0], sizeof(double), (size_t)nxny, fpd);
			fclose(fpd);
		}

		if (INDEX != 0 && INDEX % NOISEADDER == 0 && SNA == 1)
		{
			fprintf(fptmp, "Step %d: Noise added\n", INDEX);
			for (i1 = 0; i1 < n_x; ++i1)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					J = i2 + n_y * i1;
					c[J] = creal(comp[J]);
				}
			}
			add_noise(n_x, n_y, c, noise_str);
			for (i1 = 0; i1 < n_x; ++i1)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					J = i2 + n_y * i1;
					comp[J] = c[J] + _Complex_I * 0.0;
				}
			}
		}
	}
	fclose(fptmp);

	/* Free the variables and destroy the plans */

	fftw_free(g);
	fftw_free(eps_star11);
	fftw_free(eps_star12);
	fftw_free(eps_star22);
	fftw_free(u1_old);
	fftw_free(u2_old);
	fftw_free(u1_new);
	fftw_free(u2_new);
	fftw_free(mu_el);
	fftw_free(eta);

	fftw_destroy_plan(planF);
	fftw_destroy_plan(planB);

	free(c);
	free(kx);
	free(ky);
	free_dmatrix(Del_sigma_T, 1, 2, 1, 2);
	free_dmatrix(E, 1, 2, 1, 2);
	free_dmatrix(strain, 1, 2, 1, 2);
	free_dvector(Omega11, 0, nxny);
	free_dvector(Omega12, 0, nxny);
	free_dvector(Omega21, 0, nxny);
	free_dvector(Omega22, 0, nxny);
}

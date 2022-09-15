#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <complex.h>
#include "fftw3.h"

#include "../headers/functions.h"
#include "../headers/nrutil.h"

int main(void)
{

	FILE *fpr, *fpw, *fp;

	int fftw_init_threads(void);

	time_t timer;

	int n_x, n_y;
	double delta_x, delta_y;

	double kappac, kappae, A, B, P;
	double M, L;

	double eps0_11c, eps0_22c, eps0_12c;
	double eps0_11e, eps0_22e, eps0_12e;
	double c11m, c12m, c44m;
	double c11p, c12p, c44p;

	double **sig_app;

	int n_alpha, n_beta;
	int factor_alpha, factor_beta;
	int alpha_interptype, beta_interptype;
	double *alpha, *beta;
	double *alpha_prime, *beta_prime;
	double steepness_factor;

	int time_steps1, time_steps2, time_steps3;
	double delta_t1, delta_t2, delta_t3;
	double t0, t1, t2;

	int STEPS, MAXITR;
	double MAXERR;

	double ave_comp;
	fftw_complex *comp;
	fftw_complex *eta;
	double *c;
	double *e;

	double Ceff[3][3][3][3], DeltaC[3][3][3][3];
	double **sigma_Tc, **epsilon_Tc;
	double **sigma_Te, **epsilon_Te;

	int HOWTO;
	int SMOOTHEN;
	double delt;
	int tsteps;
	int INDICATOR;
	int NOISEADDER, SNA;
	double noise_str;
	int plan_indicator;
	int prob_type;

	int i, j;
	int J;
	size_t tmp;

	/* Some system calls - To get rid of old results */

	system("rm -rf ../postprocess/*.ps");
	system("rm -rf ../postprocess/*.data");
	system("rm -rf ../output/alpha*");
	system("rm -rf ../output/beta*");
	system("rm -rf ../output/comp_info");
	system("rm -rf ../output/itnos_err");
	system("rm -rf ../output/noise_info");
	system("rm -rf ../output/README");
	system("rm -rf ../output/version_name");
	system("rm -rf ../output/unconverged/*");

	printf("removed old files \n");

	/* To continue or to start anew - that is the question */

	if ((fpr = fopen("../input/howto", "r")) == NULL)
	{
		printf("Unable to open ../input/howto. Exiting\n");
		exit(0);
	}
	else
	{
		fpr = fopen("../input/howto", "r");
	}
	fscanf(fpr, "%d", &HOWTO);
	fclose(fpr);

	/* If we want to continue, let us read the relevant data and write
	 * them to the input files */

	sig_app = dmatrix(1, 2, 1, 2);
	if (HOWTO == 1)
	{

		if ((fpr = fopen("../input/continue/inputvalues", "r")) == NULL)
		{
			printf("Unable to open ../input/continue/inputvalues.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/continue/inputvalues", "r");
		}
		fscanf(fpr, "%d", &n_x);
		fscanf(fpr, "%d", &n_y);
		fscanf(fpr, "%le", &delta_x);
		fscanf(fpr, "%le", &delta_y);
		fscanf(fpr, "%le", &kappac);
		fscanf(fpr, "%le", &kappae);
		fscanf(fpr, "%le", &A);
		fscanf(fpr, "%le", &B);
		fscanf(fpr, "%le", &P);
		fscanf(fpr, "%le", &M);
		fscanf(fpr, "%le", &L);
		fscanf(fpr, "%d", &STEPS);
		fscanf(fpr, "%le", &eps0_11c);
		fscanf(fpr, "%le", &eps0_22c);
		fscanf(fpr, "%le", &eps0_12c);
		fscanf(fpr, "%le", &eps0_11e);
		fscanf(fpr, "%le", &eps0_22e);
		fscanf(fpr, "%le", &eps0_12e);
		fscanf(fpr, "%le", &c11m);
		fscanf(fpr, "%le", &c12m);
		fscanf(fpr, "%le", &c44m);
		fscanf(fpr, "%le", &c11p);
		fscanf(fpr, "%le", &c12p);
		fscanf(fpr, "%le", &c44p);
		fscanf(fpr, "%d", &n_alpha);
		fscanf(fpr, "%d", &n_beta);
		fscanf(fpr, "%d", &alpha_interptype);
		fscanf(fpr, "%d", &beta_interptype);
		fscanf(fpr, "%le", &steepness_factor);
		fscanf(fpr, "%le", &sig_app[1][1]);
		fscanf(fpr, "%le", &sig_app[2][2]);
		fscanf(fpr, "%le", &sig_app[1][2]);
		fscanf(fpr, "%d", &MAXITR);
		fscanf(fpr, "%le", &MAXERR);
		fscanf(fpr, "%le", &ave_comp);
		fscanf(fpr, "%d", &prob_type);
		fscanf(fpr, "%d", &INDICATOR);
		fscanf(fpr, "%d", &SNA);
		fscanf(fpr, "%le", &noise_str);
		fscanf(fpr, "%d", &NOISEADDER);
		fclose(fpr);

		if (alpha_interptype == 2 || alpha_interptype == 3 || alpha_interptype == 5)
		{
			printf("Please check the alpha interpolation type.\n");
			printf("In this code, they might not be relevant\n");
		}
		if (beta_interptype == 2 || beta_interptype == 3 || beta_interptype == 5)
		{
			printf("Please check the alpha interpolation type.\n");
			printf("In this code, they might not be relevant\n");
		}

		if ((fpr = fopen("../input/continue/FFTW_option", "r")) == NULL)
		{
			printf("Unable to open ../input/continue/FFTW_option.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/continue/FFTW_option", "r");
		}
		fscanf(fpr, "%d", &plan_indicator);
		fclose(fpr);

		if ((fpr = fopen("../input/continue/time_info", "r")) == NULL)
		{
			printf("Unable to open ../input/continue/time_info.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/continue/time_info", "r");
		}
		fscanf(fpr, "%le", &t0);
		fscanf(fpr, "%le%d", &delta_t1, &time_steps1);
		fscanf(fpr, "%le%d", &delta_t2, &time_steps2);
		fscanf(fpr, "%le%d", &delta_t3, &time_steps3);
		fclose(fpr);
	}

	else if (HOWTO == 0)
	{

		system("rm -rf ../output/data/*");

		if ((fpr = fopen("../input/system_constants", "r")) == NULL)
		{
			printf("Unable to open ../input/system_constants. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/system_constants", "r");
		}
		fscanf(fpr, "%d%d%le%le", &n_x, &n_y, &delta_x, &delta_y);
		fscanf(fpr, "%le%le", &kappac, &kappae);
		fscanf(fpr, "%le%le%le", &A, &B, &P);
		fscanf(fpr, "%le", &t0);
		fscanf(fpr, "%d%le", &time_steps1, &delta_t1);
		fscanf(fpr, "%d%le", &time_steps2, &delta_t2);
		fscanf(fpr, "%d%le", &time_steps3, &delta_t3);
		fscanf(fpr, "%d", &STEPS);
		fclose(fpr);

		if ((fpr = fopen("../input/mobilities", "r")) == NULL)
		{
			printf("Unable to open ../input/mobilities. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/mobilities", "r");
		}
		fscanf(fpr, "%le%le", &M, &L);
		fclose(fpr);

		if ((fpr = fopen("../input/elastic_constants", "r")) == NULL)
		{
			printf("Unable to open ../input/elastic_constants.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/elastic_constants", "r");
		}
		fscanf(fpr, "%le%le%le", &eps0_11c, &eps0_22c, &eps0_12c);
		fscanf(fpr, "%le%le%le", &eps0_11e, &eps0_22e, &eps0_12e);
		fscanf(fpr, "%le%le%le", &c11m, &c12m, &c44m);
		fscanf(fpr, "%le%le%le", &c11p, &c12p, &c44p);
		fclose(fpr);

		sig_app[2][1] = sig_app[1][2];

		if ((fpr = fopen("../input/applied_stress", "r")) == NULL)
		{
			printf("Unable to open ../input/applied_stress.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/applied_stress", "r");
		}
		fscanf(fpr, "%le%le%le", &sig_app[1][1], &sig_app[2][2], &sig_app[1][2]);
		fclose(fpr);

		if ((fpr = fopen("../input/interpolator_constants", "r")) == NULL)
		{
			printf("Unable to open\n");
			printf("../input/interpolator_constants. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/interpolator_constants", "r");
		}
		fscanf(fpr, "%d%d%d%d", &n_alpha, &n_beta, &alpha_interptype, &beta_interptype);
		fscanf(fpr, "%le", &steepness_factor);
		fclose(fpr);

		if ((fpr = fopen("../input/magic_numbers", "r")) == NULL)
		{
			printf("Unable to open ../input/magic_numbers.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/magic_numbers", "r");
		}
		fscanf(fpr, "%d%le", &MAXITR, &MAXERR);
		fclose(fpr);

		/* Sustained noise related information */

		if ((fpr = fopen("../input/sustained_noise", "r")) == NULL)
		{
			printf("Unable to open ../input/sustained_noise. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/sustained_noise", "r");
		}
		fscanf(fpr, "%d%d%le", &SNA, &NOISEADDER, &noise_str);
		fclose(fpr);

		/* Read (and write a copy of) the initial composition profile */

		if ((fpr = fopen("../input/ave_comp", "r")) == NULL)
		{
			printf("Unable to open ../input/ave_comp. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/ave_comp", "r");
		}
		tmp = fread(&ave_comp, sizeof(double), (size_t)1, fpr);
		fclose(fpr);

		/* Decide the FFTW plan flag option */

		if ((fpr = fopen("../input/FFTW_option", "r")) == NULL)
		{
			printf("Unable to open ../input/FFTW_option. Exiting");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/FFTW_option", "r");
		}
		fscanf(fpr, "%d", &plan_indicator);
		fclose(fpr);

		/* Decide if it is constrained or unconstrained evolution */

		if ((fpr = fopen("../input/evolution_type", "r")) == NULL)
		{
			printf("Unable to open ../input/evolution_type. Exiting");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/evolution_type", "r");
		}
		fscanf(fpr, "%d", &INDICATOR);
		fclose(fpr);

		/* Decide the problem type */

		if ((fpr = fopen("../input/problem_type", "r")) == NULL)
		{
			printf("Unable to open ../input/problem_type. Exiting");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/problem_type", "r");
		}
		fscanf(fpr, "%d", &prob_type);
		fclose(fpr);

	} /* End of HOWTO-else loop */

	printf("Read the input data \n");

	c = dvector(1, n_x * n_y);
	comp = fftw_malloc(n_x * n_y * sizeof(fftw_complex));

	e = dvector(1, n_x * n_y);
	eta = fftw_malloc(n_x * n_y * sizeof(fftw_complex));

	if (HOWTO == 1)
	{
		if ((fpr = fopen("../input/continue/init_compo_prof", "r")) == NULL)
		{
			printf("Unable to open ../input/continue/init_compo_prof. Exiting\n");
			exit(0);
		}
		else
		{

			fpr = fopen("../input/continue/init_compo_prof", "r");
		}
		tmp = fread(&c[1], sizeof(double), (size_t)n_x * n_y, fpr);
		fclose(fpr);
	}
	else if (HOWTO == 0)
	{
		/** Read compsotion **/

		if ((fpr = fopen("../input/init_compo_prof", "r")) == NULL)
		{
			printf("Unable to open ../input/init_compo_prof. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/init_compo_prof", "r");
		}
		tmp = fread(&c[1], sizeof(double), (size_t)n_x * n_y, fpr);
		fclose(fpr);

		/** Read eta **/

		if ((fpr = fopen("../input/init_eta_prof", "r")) == NULL)
		{
			printf("Unable to open ../input/init_eta_prof. Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../input/init_eta_prof", "r");
		}
		tmp = fread(&e[1], sizeof(double), (size_t)n_x * n_y, fpr);
		fclose(fpr);
	}
	else
	{
		printf("It is not clear whether I should start anew or continue\n");
		printf("an older calculation. Check the input file input/howto\n");
		printf("I am exiting.\n");
		exit(0);
	}
	for (i = 0; i < n_x; ++i)
	{
		for (j = 0; j < n_y; ++j)
		{
			J = j + n_y * i;
			__real__ comp[J] = c[J + 1];
			__imag__ comp[J] = 0.0;

			__real__ eta[J] = e[J + 1];
			__imag__ eta[J] = 0.0;
		}
	}

	if ((fpw = fopen("../output/comp_info", "w")) == NULL)
	{
		printf("Unable to open ../output/comp_info.");
		printf("Exiting\n");
		exit(0);
	}
	else
	{
		fpw = fopen("../output/comp_info", "w");
	}
	fprintf(fpw, "Average composition: %le\n", ave_comp);
	fprintf(fpw, "Composition profile:\n");
	for (i = 0; i < n_x; ++i)
	{
		for (j = 0; j < n_y; ++j)
		{
			fprintf(fpw, "%d %d %le\n", i, j, __real__ comp[j + n_y * i]);
		}
	}
	fclose(fpw);

	if ((fpw = fopen("../output/eta_info", "w")) == NULL)
	{
		printf("Unable to open ../output/eta_info.");
		printf("Exiting\n");
		exit(0);
	}
	else
	{
		fpw = fopen("../output/eta_info", "w");
	}
	fprintf(fpw, "Eta profile:\n");
	for (i = 0; i < n_x; ++i)
	{
		for (j = 0; j < n_y; ++j)
		{
			fprintf(fpw, "%d %d %le\n", i, j, __real__ eta[j + n_y * i]);
		}
	}
	fclose(fpw);

	printf("Written the comp and eta profiles \n");

	/* Generate the elasticity matrices */

	generate_elasticity_tensor(c11m, c12m, c44m, c11p, c12p, c44p, Ceff, DeltaC,
							   alpha_interptype, prob_type, ave_comp);

	printf("Generated elasticity tensor\n");

	/* Calculate the eigen (T as in transformational) stress and strain */

	sigma_Tc = dmatrix(1, 2, 1, 2);
	epsilon_Tc = dmatrix(1, 2, 1, 2);
	sigma_Te = dmatrix(1, 2, 1, 2);
	epsilon_Te = dmatrix(1, 2, 1, 2);
	calculate_sigma_T(Ceff, eps0_11c, eps0_22c, eps0_12c, sigma_Tc, epsilon_Tc);
	calculate_sigma_T(Ceff, eps0_11e, eps0_22e, eps0_12e, sigma_Te, epsilon_Te);

	printf("Calculated eigen strain and stress\n");

	system("rm -rf ../output/SG_Analysis");
	if (sig_app[1][1] != 0.0 || sig_app[2][2] != 0.0)
	{
		calculate_aleph(c11m, c12m, c44m, c11p, c12p, c44p, sig_app, epsilon_Tc);
	}

	/* Calculate the time translations and write the postprocessor files */

	t1 = time_steps1 + t0;
	t2 = time_steps2 + t1;
	system("rm -rf ../output/postprocessor/ps/src/psmaker.c");
	system("rm -rf ../output/postprocessor/data/src/datamaker.c");
	write_ppfiles(n_x, n_y, delta_t1, time_steps1, t0, delta_t2, time_steps2, t1,
				  delta_t3, time_steps3, t2, STEPS);

	printf("Calculated time translations\n");

	/* Keep a copy of input values, for use in case the program is
	 * terminated unexpectedly */

	system("rm -rf ../output/cotinue/*");
	if ((fpw = fopen("../output/continue/README", "w")) == NULL)
	{
		printf("Unable to open ../output/continue/README.");
		printf("Exiting\n");
		exit(0);
	}
	else
	{
		fpw = fopen("../output/continue/README", "w");
	}
	fprintf(fpw, "%d\n", n_x);
	fprintf(fpw, "%d\n", n_y);
	fprintf(fpw, "%le\n", delta_x);
	fprintf(fpw, "%le\n", delta_y);
	fprintf(fpw, "%le\n", kappac);
	fprintf(fpw, "%le\n", kappae);
	fprintf(fpw, "%le\n", A);
	fprintf(fpw, "%le\n", B);
	fprintf(fpw, "%le\n", P);
	fprintf(fpw, "%le\n", M);
	fprintf(fpw, "%le\n", L);
	fprintf(fpw, "%d\n", STEPS);
	fprintf(fpw, "%le\n", eps0_11c);
	fprintf(fpw, "%le\n", eps0_22c);
	fprintf(fpw, "%le\n", eps0_12c);
	fprintf(fpw, "%le\n", eps0_11e);
	fprintf(fpw, "%le\n", eps0_22e);
	fprintf(fpw, "%le\n", eps0_12e);
	fprintf(fpw, "%le\n", c11m);
	fprintf(fpw, "%le\n", c12m);
	fprintf(fpw, "%le\n", c44m);
	fprintf(fpw, "%le\n", c11p);
	fprintf(fpw, "%le\n", c12p);
	fprintf(fpw, "%le\n", c44p);
	fprintf(fpw, "%d\n", n_alpha);
	fprintf(fpw, "%d\n", n_beta);
	fprintf(fpw, "%d\n", alpha_interptype);
	fprintf(fpw, "%d\n", beta_interptype);
	fprintf(fpw, "%le\n", steepness_factor);
	fprintf(fpw, "%le\n", sig_app[1][1]);
	fprintf(fpw, "%le\n", sig_app[2][2]);
	fprintf(fpw, "%le\n", sig_app[1][2]);
	fprintf(fpw, "%d\n", MAXITR);
	fprintf(fpw, "%le\n", MAXERR);
	fprintf(fpw, "%le\n", ave_comp);
	fprintf(fpw, "%d\n", prob_type);
	fprintf(fpw, "%d\n", INDICATOR);
	fprintf(fpw, "%d\n", SNA);
	fprintf(fpw, "%le\n", noise_str);
	fprintf(fpw, "%d\n", NOISEADDER);
	fclose(fpw);

	factor_alpha = (int)n_alpha;
	factor_beta = (int)n_beta;
	alpha = dvector(-factor_alpha, n_alpha + factor_alpha);
	beta = dvector(-factor_beta, n_alpha + factor_beta);
	alpha_prime = dvector(-factor_alpha, n_alpha + factor_alpha);
	beta_prime = dvector(-factor_beta, n_alpha + factor_beta);

	/* Open the README file - To write the input values used */

	printf("Opening README file\n");

	if ((fpw = fopen("../output/README", "w")) == NULL)
	{
		printf("Unable to open ../output/README.");
		printf("Exiting\n");
		exit(0);
	}
	else
	{
		fpw = fopen("../output/README", "w");
	}

	timer = time(NULL);

	fprintf(fpw, "Data in this directory is generated on");
	fprintf(fpw, " %s", asctime(localtime(&timer)));
	fprintf(fpw, "See the file version_name in this directory\n");
	fprintf(fpw, "to find the version of the code used\n");
	system("pwd > ../output/version_name");
	if (HOWTO == 0)
		fprintf(fpw, "Starting the calculations afresh\n");
	else if (HOWTO == 1)
		fprintf(fpw, "Continuation of an old calculation.\n");

	/* To check if the initial profile needs to be smoothened */

	if ((fpr = fopen("../input/smoothen", "r")) == NULL)
	{
		printf("Unable to open ../input/smoothen. Exiting\n");
		exit(0);
	}
	else
	{
		fpr = fopen("../input/smoothen", "r");
	}
	fscanf(fpr, "%d", &SMOOTHEN);
	fscanf(fpr, "%le%d", &delt, &tsteps);
	fclose(fpr);

	/* If needed, smoothen and write the smoothened profile to output */

	if (SMOOTHEN == 1)
	{

		evolve_diffusional(n_x, n_y, delta_x, delta_y, kappac, A, delt, tsteps, 0,
						   comp, ave_comp, tsteps + 100, plan_indicator, 0, tsteps + 100, 0.0);

		system("rm -rf ../output/data/*");

		if ((fpr = fopen("../output/smoothened_comp_info", "w")) == NULL)
		{
			printf("Unable to open ../output/smoothened_comp_info.");
			printf("Exiting\n");
			exit(0);
		}
		else
		{
			fpr = fopen("../output/smoothened_comp_info", "w");
		}
		fprintf(fpr, "Average composition: %le\n", ave_comp);
		fprintf(fpr, "Composition profile:\n");
		for (i = 0; i < n_x; ++i)
		{
			for (j = 0; j < n_y; ++j)
			{
				fprintf(fpr, "%d %d %le\n", i, j, __real__(comp[j + n_y * i]));
				if (__real__ comp[j + n_y * i] >= 1.5)
					printf("1 %d %d %le\n", i, j, __real__ comp[j + n_y * i]);
			}
		}
		fclose(fpr);

		/* Put a note in the README file */

		fprintf(fpw, "The initial profile was smoothened:\n");
		fprintf(fpw, "for %d time steps of delta_t %le\n", tsteps, delt);
	}

	/* Evolve the microstructure */

	if (prob_type == 1)
	{

		fprintf(fpw, "Diffusional evolution without elasticity\n");
		fprintf(fpw, "n_x: %d n_y: %d\n", n_x, n_y);
		fprintf(fpw, "delta_x: %le delta_y: %le\n", delta_x, delta_y);
		fprintf(fpw, "kappa: %le A: %le\n", kappac, A);
		fprintf(fpw, "t0: %le\n", t0);
		fprintf(fpw, "time_steps1: %d delta_t1: %le\n", time_steps1, delta_t1);
		fprintf(fpw, "time_steps2: %d delta_t2: %le\n", time_steps2, delta_t2);
		fprintf(fpw, "time_steps3: %d delta_t3: %le\n", time_steps3, delta_t3);
		fprintf(fpw, "STEPS: %d\n", STEPS);
		fprintf(fpw, "Alloy composition: %le\n", ave_comp);
		if (SNA == 0)
		{
			fprintf(fpw, "No sustained noise.\n");
		}
		else if (SNA == 1)
		{
			fprintf(fpw, "Sustained noise of strength %le added:", noise_str);
			fprintf(fpw, "Once in %d steps\n", NOISEADDER);
		}
		fclose(fpw);

		evolve_diffusional(n_x, n_y, delta_x, delta_y, kappac, A, delta_t1, time_steps1, t0,
						   comp, ave_comp, STEPS, plan_indicator, SNA, NOISEADDER, noise_str);
		evolve_diffusional(n_x, n_y, delta_x, delta_y, kappac, A, delta_t2, time_steps2, t1,
						   comp, ave_comp, STEPS, plan_indicator, SNA, NOISEADDER, noise_str);
		evolve_diffusional(n_x, n_y, delta_x, delta_y, kappac, A, delta_t3, time_steps3, t2,
						   comp, ave_comp, STEPS, plan_indicator, SNA, NOISEADDER, noise_str);
	}

	if (prob_type == 2)
	{

		/* Read (and write a copy of) the interpolation functions */

		calculate_beta(n_beta, beta, beta_interptype, steepness_factor, ave_comp);
		calculate_beta_prime(n_beta, beta_prime, beta_interptype, steepness_factor);

		if ((fp = fopen("../output/beta", "w")) == NULL)
		{
			printf("Unable to open ../output/beta. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/beta", "w");
		}

		for (i = -factor_beta; i < n_beta + factor_beta + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, beta[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/beta_prime", "w")) == NULL)
		{
			printf("Unable to open ../output/beta_prime. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/beta_prime", "w");
		}

		for (i = -factor_beta; i < n_beta + factor_beta + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, beta_prime[i]);
		}

		fclose(fp);

		fprintf(fpw, "Evolution with homogeneous elasticity\n");

		if (INDICATOR == 0)
		{
			fprintf(fpw, "Constrained evolution\n");
		}
		else if (INDICATOR == 1)
		{
			fprintf(fpw, "Unconstrained evolution\n");
		}

		fprintf(fpw, "n_x: %d n_y: %d\n", n_x, n_y);
		fprintf(fpw, "delta_x: %le delta_y: %le\n", delta_x, delta_y);
		fprintf(fpw, "kappa: %le A: %le\n", kappac, A);
		fprintf(fpw, "t0: %le\n", t0);
		fprintf(fpw, "time_steps1: %d delta_t1: %le\n", time_steps1, delta_t1);
		fprintf(fpw, "time_steps2: %d delta_t2: %le\n", time_steps2, delta_t2);
		fprintf(fpw, "time_steps3: %d delta_t3: %le\n", time_steps3, delta_t3);
		fprintf(fpw, "STEPS: %d\n", STEPS);
		fprintf(fpw, "eps0_11: %le\t", eps0_11c);
		fprintf(fpw, "eps0_22: %le\t", eps0_22c);
		fprintf(fpw, "eps0_12: %le\n", eps0_12c);
		fprintf(fpw, "c11: %le\t", c11m);
		fprintf(fpw, "c12: %le\t", c12m);
		fprintf(fpw, "c44: %le\n", c44m);
		fprintf(fpw, "n_beta: %d\n", n_beta);
		fprintf(fpw, "beta_interptype: %d:", beta_interptype);
		if (beta_interptype == 1)
			fprintf(fpw, "This means it is a tanh interpolation\n");
		else if (beta_interptype == 2)
			fprintf(fpw, "This means it is a comp interpolation\n");
		else if (beta_interptype == 3)
			fprintf(fpw, "This means it is a del_comp interpolation\n");
		else if (beta_interptype == 4)
			fprintf(fpw, "This means it is an interpolation by Wang's function\n");
		fprintf(fpw, "sig_app[1][1]: %le\t", sig_app[1][1]);
		fprintf(fpw, "sig_app[2][2]: %le\n", sig_app[2][2]);
		fprintf(fpw, "sig_app[1][2]: %le\t", sig_app[1][2]);

		fprintf(fpw, "Alloy composition: %le\n", ave_comp);

		if (SNA == 0)
		{
			fprintf(fpw, "No sustained noise.\n");
		}
		else if (SNA == 1)
		{
			fprintf(fpw, "Sustained noise of strength %le added:", noise_str);
			fprintf(fpw, "Once in %d steps\n", NOISEADDER);
		}
		fclose(fpw);

		if (INDICATOR == 0)
		{
			evolve_hom_constrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t1, time_steps1,
								   t0, comp, Ceff, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS, n_beta, beta,
								   beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_hom_constrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t2, time_steps2,
								   t1, comp, Ceff, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS, n_beta, beta,
								   beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_hom_constrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t3, time_steps3,
								   t2, comp, Ceff, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS, n_beta, beta,
								   beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
		}
		else
		{
			evolve_hom_unconstrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t1, time_steps1,
									 t0, comp, Ceff, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									 n_beta, beta, beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_hom_unconstrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t2, time_steps2,
									 t1, comp, Ceff, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									 n_beta, beta, beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_hom_unconstrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t3, time_steps3,
									 t2, comp, Ceff, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									 n_beta, beta, beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
		}
	}

	if (prob_type == 3)
	{

		printf("Reached prob_type 3\n");

		/* Read (and write a copy of) the interpolation functions */

		calculate_alpha(n_alpha, alpha, alpha_interptype, steepness_factor, ave_comp);
		calculate_beta(n_beta, beta, beta_interptype, steepness_factor, ave_comp);
		calculate_alpha_prime(n_alpha, alpha_prime, alpha_interptype, steepness_factor);
		calculate_beta_prime(n_beta, beta_prime, beta_interptype, steepness_factor);

		printf("Finished Interpolation function calculations\n");

		if ((fp = fopen("../output/alpha", "w")) == NULL)
		{
			printf("Unable to open ../output/alpha. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/alpha", "w");
		}

		for (i = -factor_alpha; i < n_alpha + factor_alpha + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, alpha[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/beta", "w")) == NULL)
		{
			printf("Unable to open ../output/beta. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/beta", "w");
		}

		for (i = -factor_beta; i < n_beta + factor_beta + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, beta[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/alpha_prime", "w")) == NULL)
		{
			printf("Unable to open ../output/alpha_prime. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/alpha_prime", "w");
		}

		for (i = -factor_alpha; i < n_alpha + factor_alpha + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, alpha_prime[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/beta_prime", "w")) == NULL)
		{
			printf("Unable to open ../output/beta_prime. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/beta_prime", "w");
		}

		for (i = -factor_beta; i < n_beta + factor_beta + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, beta_prime[i]);
		}

		fclose(fp);

		fprintf(fpw, "Evolution with inhomogeneous elasticity\n");
		if (INDICATOR == 0)
		{
			fprintf(fpw, "Constrained evolution\n");
		}
		else if (INDICATOR == 1)
		{
			fprintf(fpw, "Unconstrained evolution\n");
		}

		fprintf(fpw, "n_x: %d n_y: %d\n", n_x, n_y);
		fprintf(fpw, "delta_x: %le delta_y: %le\n", delta_x, delta_y);
		fprintf(fpw, "kappac: %le kappae: %le\n", kappac, kappae);
		fprintf(fpw, "A: %le B: %le P: %le\n", A, B, P);
		fprintf(fpw, "t0: %le\n", t0);
		fprintf(fpw, "time_steps1: %d delta_t1: %le\n", time_steps1, delta_t1);
		fprintf(fpw, "time_steps2: %d delta_t2: %le\n", time_steps2, delta_t2);
		fprintf(fpw, "time_steps3: %d delta_t3: %le\n", time_steps3, delta_t3);
		fprintf(fpw, "STEPS: %d\n", STEPS);
		fprintf(fpw, "eps0_11c: %le\t", eps0_11c);
		fprintf(fpw, "eps0_22c: %le\t", eps0_22c);
		fprintf(fpw, "eps0_12c: %le\n", eps0_12c);
		fprintf(fpw, "eps0_11e: %le\t", eps0_11e);
		fprintf(fpw, "eps0_22e: %le\t", eps0_22e);
		fprintf(fpw, "eps0_12e: %le\n", eps0_12e);
		fprintf(fpw, "c11m: %le\t", c11m);
		fprintf(fpw, "c12m: %le\t", c12m);
		fprintf(fpw, "c44m: %le\n", c44m);
		fprintf(fpw, "c11p: %le\t", c11p);
		fprintf(fpw, "c12p: %le\t", c12p);
		fprintf(fpw, "c44p: %le\n", c44p);
		fprintf(fpw, "n_alpha: %d\t", n_alpha);
		fprintf(fpw, "n_beta: %d\n", n_beta);
		fprintf(fpw, "alpha_interptype: %d:", alpha_interptype);
		if (alpha_interptype == 1)
			fprintf(fpw, "This means it is a tanh interpolation\n");
		else if (alpha_interptype == 2)
			fprintf(fpw, "This means it is a comp interpolation\n");
		else if (alpha_interptype == 3)
			fprintf(fpw, "This means it is a del_comp interpolation\n");
		else if (alpha_interptype == 4)
			fprintf(fpw, "This means it is an interpolation by Wang's function\n");
		else if (alpha_interptype == 5)
			fprintf(fpw, "This means it is a comp-1 interpolation\n");
		fprintf(fpw, "beta_interptype: %d:", beta_interptype);
		if (beta_interptype == 1)
			fprintf(fpw, "This means it is a tanh interpolation\n");
		else if (beta_interptype == 2)
			fprintf(fpw, "This means it is a comp interpolation\n");
		else if (beta_interptype == 3)
			fprintf(fpw, "This means it is a del_comp interpolation\n");
		else if (beta_interptype == 4)
			fprintf(fpw, "This means it is an interpolation by Wang's function\n");
		else if (beta_interptype == 5)
			fprintf(fpw, "This means it is a comp interpolation\n");
		fprintf(fpw, "sig_app[1][1]: %le\t", sig_app[1][1]);
		fprintf(fpw, "sig_app[2][2]: %le\n", sig_app[2][2]);
		fprintf(fpw, "sig_app[1][2]: %le\t", sig_app[1][2]);
		fprintf(fpw, "Maximum number of iterations allowed: %d and \n", MAXITR);
		fprintf(fpw, "maximum allowed error in displacement refinement: %le\n", MAXERR);
		fprintf(fpw, "Alloy composition: %le\n", ave_comp);

		if (SNA == 0)
		{
			fprintf(fpw, "No sustained noise.\n");
		}
		else if (SNA == 1)
		{
			fprintf(fpw, "Sustained noise of strength %le added:", noise_str);
			fprintf(fpw, "Once in %d steps\n", NOISEADDER);
		}
		fclose(fpw);
		//
		for (i = 0; i < n_x; ++i)
		{
			for (j = 0; j < n_y; ++j)
			{
				if (__real__ comp[j + n_y * i] >= 1.5)
					printf("2 %d %d %le\n", i, j, __real__ comp[j + n_y * i]);
			}
		}
		//
		if (INDICATOR == 0)
		{
			evolve_inhom_constrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t1,
									 time_steps1, t0, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									 MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime, beta_prime, plan_indicator,
									 SNA, NOISEADDER, noise_str);
			evolve_inhom_constrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t2,
									 time_steps2, t1, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									 MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime, beta_prime, plan_indicator,
									 SNA, NOISEADDER, noise_str);
			evolve_inhom_constrained(n_x, n_y, delta_x, delta_y, kappac, A, delta_t3,
									 time_steps3, t2, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									 MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime, beta_prime, plan_indicator,
									 SNA, NOISEADDER, noise_str);
		}
		else
		{
			/*evolve_inhom_unconstrained_ceta(n_x,n_y,delta_x,delta_y,kappac,kappae,
			A,B,P,M,L,delta_t1,time_steps1,t0,comp,eta,Ceff,DeltaC,sigma_Tc,epsilon_Tc,
			sigma_Te,epsilon_Te,ave_comp,sig_app,STEPS,MAXITR,MAXERR,n_alpha,n_beta,
			alpha,beta,alpha_prime,beta_prime,plan_indicator,SNA,NOISEADDER,noise_str);
			evolve_inhom_unconstrained_ceta(n_x,n_y,delta_x,delta_y,kappac,kappae,
			A,B,P,M,L,delta_t2,time_steps2,t1,comp,eta,Ceff,DeltaC,sigma_Tc,epsilon_Tc,
			sigma_Te,epsilon_Te,ave_comp,sig_app,STEPS,MAXITR,MAXERR,n_alpha,n_beta,
			alpha,beta,alpha_prime,beta_prime,plan_indicator,SNA,NOISEADDER,noise_str);
			evolve_inhom_unconstrained_ceta(n_x,n_y,delta_x,delta_y,kappac,kappae,
			A,B,P,M,L,delta_t3,time_steps3,t2,comp,eta,Ceff,DeltaC,sigma_Tc,epsilon_Tc,
			sigma_Te,epsilon_Te,ave_comp,sig_app,STEPS,MAXITR,MAXERR,n_alpha,n_beta,
			alpha,beta,alpha_prime,beta_prime,plan_indicator,SNA,NOISEADDER,noise_str);*/

			evolve_inhom_unconstrained(n_x, n_y, delta_x, delta_y, kappac,
									   A, delta_t1, time_steps1, t0, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc,
									   ave_comp, sig_app, STEPS, MAXITR, MAXERR, n_alpha, n_beta,
									   alpha, beta, alpha_prime, beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_inhom_unconstrained(n_x, n_y, delta_x, delta_y, kappac,
									   A, delta_t2, time_steps2, t1, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc,
									   ave_comp, sig_app, STEPS, MAXITR, MAXERR, n_alpha, n_beta,
									   alpha, beta, alpha_prime, beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_inhom_unconstrained(n_x, n_y, delta_x, delta_y, kappac,
									   A, delta_t3, time_steps3, t2, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc,
									   ave_comp, sig_app, STEPS, MAXITR, MAXERR, n_alpha, n_beta,
									   alpha, beta, alpha_prime, beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
		}
	}
	if (prob_type == 4)
	{

		/* Read (and write a copy of) the interpolation functions */

		calculate_alpha(n_alpha, alpha, alpha_interptype, steepness_factor, ave_comp);
		calculate_beta(n_beta, beta, beta_interptype, steepness_factor, ave_comp);
		calculate_alpha_prime(n_alpha, alpha_prime, alpha_interptype, steepness_factor);
		calculate_beta_prime(n_beta, beta_prime, beta_interptype, steepness_factor);

		if ((fp = fopen("../output/alpha", "w")) == NULL)
		{
			printf("Unable to open ../output/alpha. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/alpha", "w");
		}

		for (i = -factor_alpha; i < n_alpha + factor_alpha + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, alpha[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/beta", "w")) == NULL)
		{
			printf("Unable to open ../output/beta. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/beta", "w");
		}

		for (i = -factor_beta; i < n_beta + factor_beta + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, beta[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/alpha_prime", "w")) == NULL)
		{
			printf("Unable to open ../output/alpha_prime. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/alpha_prime", "w");
		}

		for (i = -factor_alpha; i < n_alpha + factor_alpha + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, alpha_prime[i]);
		}

		fclose(fp);

		if ((fp = fopen("../output/beta_prime", "w")) == NULL)
		{
			printf("Unable to open ../output/beta_prime. Exiting.\n");
			exit(0);
		}
		else
		{
			fp = fopen("../output/beta_prime", "w");
		}

		for (i = -factor_beta; i < n_beta + factor_beta + 1; ++i)
		{
			fprintf(fp, "%d %le\n", i, beta_prime[i]);
		}

		fclose(fp);

		fprintf(fpw, "Evolution with\n");
		fprintf(fpw, "inhomogeneous elasticity and variable mobility\n");
		if (INDICATOR == 0)
		{
			fprintf(fpw, "Constrained evolution\n");
		}
		else if (INDICATOR == 1)
		{
			fprintf(fpw, "Unconstrained evolution\n");
		}

		fprintf(fpw, "n_x: %d n_y: %d\n", n_x, n_y);
		fprintf(fpw, "delta_x: %le delta_y: %le\n", delta_x, delta_y);
		fprintf(fpw, "kappa: %le A: %le\n", kappac, A);
		fprintf(fpw, "M_b: %le M_n: %le M_t: %le\n", M, L, L);
		fprintf(fpw, "t0: %le\n", t0);
		fprintf(fpw, "time_steps1: %d delta_t1: %le\n", time_steps1, delta_t1);
		fprintf(fpw, "time_steps2: %d delta_t2: %le\n", time_steps2, delta_t2);
		fprintf(fpw, "time_steps3: %d delta_t3: %le\n", time_steps3, delta_t3);
		fprintf(fpw, "STEPS: %d\n", STEPS);
		fprintf(fpw, "eps0_11: %le\t", eps0_11c);
		fprintf(fpw, "eps0_22: %le\t", eps0_22c);
		fprintf(fpw, "eps0_12: %le\n", eps0_12c);
		fprintf(fpw, "c11m: %le\t", c11m);
		fprintf(fpw, "c12m: %le\t", c12m);
		fprintf(fpw, "c44m: %le\n", c44m);
		fprintf(fpw, "c11p: %le\t", c11p);
		fprintf(fpw, "c12p: %le\t", c12p);
		fprintf(fpw, "c44p: %le\n", c44p);
		fprintf(fpw, "n_alpha: %d\t", n_alpha);
		fprintf(fpw, "n_beta: %d\n", n_beta);
		fprintf(fpw, "alpha_interptype: %d:", alpha_interptype);
		if (alpha_interptype == 1)
			fprintf(fpw, "This means it is a tanh interpolation\n");
		else if (alpha_interptype == 2)
			fprintf(fpw, "This means it is a comp interpolation\n");
		else if (alpha_interptype == 3)
			fprintf(fpw, "This means it is a del_comp interpolation\n");
		else if (alpha_interptype == 4)
			fprintf(fpw, "This means it is an interpolation by Wang's function\n");
		else if (alpha_interptype == 5)
			fprintf(fpw, "This means it is a comp-1 interpolation\n");
		fprintf(fpw, "beta_interptype: %d:", beta_interptype);
		if (beta_interptype == 1)
			fprintf(fpw, "This means it is a tanh interpolation\n");
		else if (beta_interptype == 2)
			fprintf(fpw, "This means it is a comp interpolation\n");
		else if (beta_interptype == 3)
			fprintf(fpw, "This means it is a del_comp interpolation\n");
		else if (beta_interptype == 4)
			fprintf(fpw, "This means it is an interpolation by Wang's function\n");
		else if (beta_interptype == 5)
			fprintf(fpw, "This means it is a comp interpolation\n");
		fprintf(fpw, "sig_app[1][1]: %le\t", sig_app[1][1]);
		fprintf(fpw, "sig_app[2][2]: %le\n", sig_app[2][2]);
		fprintf(fpw, "sig_app[1][2]: %le\t", sig_app[1][2]);
		fprintf(fpw, "Maximum number of iterations allowed: %d and \n", MAXITR);
		fprintf(fpw, "maximum allowed error in displacement refinement: %le\n", MAXERR);
		fprintf(fpw, "Alloy composition: %le\n", ave_comp);

		if (SNA == 0)
		{
			fprintf(fpw, "No sustained noise.\n");
		}
		else if (SNA == 1)
		{
			fprintf(fpw, "Sustained noise of strength %le added:", noise_str);
			fprintf(fpw, "Once in %d steps\n", NOISEADDER);
		}
		fclose(fpw);

		if (INDICATOR == 0)
		{
			evolve_inhom_const_varmob(n_x, n_y, delta_x, delta_y, kappac, A, M, L, L, delta_t1,
									  time_steps1, t0, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									  MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime, beta_prime, plan_indicator,
									  SNA, NOISEADDER, noise_str);
			evolve_inhom_const_varmob(n_x, n_y, delta_x, delta_y, kappac, A, M, L, L, delta_t2,
									  time_steps2, t1, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									  MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime, beta_prime, plan_indicator,
									  SNA, NOISEADDER, noise_str);
			evolve_inhom_const_varmob(n_x, n_y, delta_x, delta_y, kappac, A, M, L, L, delta_t3,
									  time_steps3, t2, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp, sig_app, STEPS,
									  MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime, beta_prime, plan_indicator,
									  SNA, NOISEADDER, noise_str);
		}
		else
		{
			evolve_inhom_unconst_varmob(n_x, n_y, delta_x, delta_y, kappac, A, M, L, L,
										delta_t1, time_steps1, t0, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp,
										sig_app, STEPS, MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime,
										beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_inhom_unconst_varmob(n_x, n_y, delta_x, delta_y, kappac, A, M, L, L,
										delta_t2, time_steps2, t1, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp,
										sig_app, STEPS, MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime,
										beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
			evolve_inhom_unconst_varmob(n_x, n_y, delta_x, delta_y, kappac, A, M, L, L,
										delta_t3, time_steps3, t2, comp, Ceff, DeltaC, sigma_Tc, epsilon_Tc, ave_comp,
										sig_app, STEPS, MAXITR, MAXERR, n_alpha, n_beta, alpha, beta, alpha_prime,
										beta_prime, plan_indicator, SNA, NOISEADDER, noise_str);
		}
	}

	/* Clean up */

	fftw_free(comp);
	free_dmatrix(sig_app, 1, 2, 1, 2);
	free_dmatrix(sigma_Tc, 1, 2, 1, 2);
	free_dmatrix(epsilon_Tc, 1, 2, 1, 2);
	free_dvector(c, 1, n_x * n_y);
	free_dvector(alpha, -factor_alpha, n_alpha + factor_alpha);
	free_dvector(beta, -factor_beta, n_beta + factor_beta);
	free_dvector(alpha_prime, -factor_alpha, n_alpha + factor_alpha);
	free_dvector(beta_prime, -factor_beta, n_beta + factor_beta);

	/* Say bye */

	printf("Program run successfully. Time to say adieu\n");

	void fftw_cleanup_threads(void);

	return 0;
}

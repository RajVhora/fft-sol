/**************************************************************
IF YOU ARE MODIFYING THIS CODE, YOU SHOULD ALSO MODIFY
calculate_beta. The following table indicates the interpolation
type
1 :  tanh(c)
2 :  c
3 :  c-ave_comp
4 :  Wang's function
5 :  c
***************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void calculate_beta_prime(int n_beta, double *beta, 
int beta_interptype, double steepness_factor){

int i;

double onebyf;
double jnk;
int factor;
double bp;

onebyf = 1.0/n_beta;
jnk = -1.0;
factor = (int) (n_beta);

if(beta_interptype == 1){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_beta+factor+1; ++i){
	  bp = 1.0/(cosh(2.0*steepness_factor*jnk-steepness_factor));
		beta[i] = 	steepness_factor*bp*bp;
			if(jnk < 0.0) beta[i] = 0.0;
			else if(jnk > 1.0) beta[i] = 0.0;
		jnk = jnk + onebyf;
	}
}
else if(beta_interptype == 2 || beta_interptype == 5){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_beta+factor+1; ++i){
		beta[i] = 	1.0;
		jnk = jnk + onebyf;
	}
}
else if(beta_interptype == 3){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_beta+factor+1; ++i){
		beta[i] = 	1.0;
		jnk = jnk + onebyf;
	}
}
else if(beta_interptype == 4){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_beta+factor+1; ++i){
		beta[i] = 	30.0*jnk*jnk*(1.0-2.0*jnk+jnk*jnk);
			if(jnk < 0.0) beta[i] = 0.0;
			else if(jnk > 1.0) beta[i] = 0.0;
		jnk = jnk + onebyf;
	}
}
else{
printf("Unable to determine beta interpolation type\n");
printf("Exiting from calculate_beta_prime.c\n");
exit(0);
}

}


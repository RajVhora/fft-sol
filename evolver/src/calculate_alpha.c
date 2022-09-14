/********************************************************************* 
IF YOU ARE MODIFYING THIS CODE, YOU SHOULD ALSO MODIFY
calculate_alpha_prime. The following table indicates the interpolation
type
1 :  tanh(c)
2 :  c
3 :  c-ave_comp
4 :  Wang's function
5 :  c-1
**********************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void calculate_alpha(int n_alpha, double *alpha, int alpha_interptype,
double steepness_factor, double ave_comp){

int i;

double onebyf;
double jnk;
int factor;

onebyf = 1.0/n_alpha;
jnk = -1.0;
factor = (int) (n_alpha);

if(alpha_interptype == 1){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_alpha+factor+1; ++i){
		alpha[i] = 	0.5*(tanh(2.0*steepness_factor*jnk-steepness_factor));
			if(jnk < 0.0) alpha[i] = -0.5;
			else if(jnk > 1.0) alpha[i] = 0.5;
		jnk = jnk + onebyf;
	}
}
else if(alpha_interptype == 2){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_alpha+factor+1; ++i){
		alpha[i] = 	jnk;
			if(jnk < 0.0) alpha[i] = 0.0;
			else if(jnk > 1.0) alpha[i] = 1.0;
		jnk = jnk + onebyf;
	}
}
else if(alpha_interptype == 3){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_alpha+factor+1; ++i){
		alpha[i] = 	jnk - ave_comp;
		jnk = jnk + onebyf;
	}
}
else if(alpha_interptype == 4){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_alpha+factor+1; ++i){
		alpha[i] = 	jnk*jnk*jnk*(10.0-15.0*jnk+6.0*jnk*jnk) - 0.5;
			if(jnk < 0.0) alpha[i] = -0.5;
			else if(jnk > 1.0) alpha[i] = 0.5;
		jnk = jnk + onebyf;
	}
}
else if(alpha_interptype == 5){
/* Before editing this section, see the beginning of this file */
	for(i=-factor; i < n_alpha+factor+1; ++i){
		alpha[i] = 	jnk - 1.0;
			if(jnk < 0.0) alpha[i] = -1.0;
			else if(jnk > 1.0) alpha[i] = 0.0;
		jnk = jnk + onebyf;
	}
}
else{
printf("Unable to determine the alpha interpolation type.\n");
printf("Exiting from calculate_alpha.c\n");
exit(0);
}

}


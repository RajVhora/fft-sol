#include<math.h>

void elliptic(int n_x, int n_y, double R1, double R2,
double comp_precipitate, double comp_matrix, double *comp){

int i1,i2;
double a,b,c;
double tmp1,tmp2;
double n_x_by2,n_y_by2;

R1 = R1*R1;
R2 = R2*R2;
a = R1;
b = R2;
n_x_by2 = (int) (n_x/2) - 1;
n_y_by2 = (int) (n_y/2) - 1;
for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_y; ++i2){
	tmp1 = 1.0*i1 - n_x/2.0;
	tmp1 = tmp1*tmp1/a;
	tmp2 = 1.0*i2 - n_y/2.0;
	tmp2 = tmp2*tmp2/b;
	c = tmp1 + tmp2;
	if( c <= 1.0 )
		comp[i2+n_y*i1] = comp_precipitate;
	else
		comp[i2+n_y*i1] = comp_matrix;
}}

}

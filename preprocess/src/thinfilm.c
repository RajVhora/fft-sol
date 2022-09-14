#include<math.h>

void thinfilm(int n_x, int n_y, int h, double comp_precipitate,
double comp_matrix, double *comp){

int i1,i2;
int n_y_by2;
int halfh;

n_y_by2 = (int) (n_y/2) - 1;
halfh = (int) (h/2);

for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_y; ++i2){
	if(i2 >= n_y_by2-halfh && i2 <= n_y_by2+halfh)
		comp[i2+n_y*i1] = comp_precipitate;
	else
		comp[i2+n_y*i1] = comp_matrix;
}}

}

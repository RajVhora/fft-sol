
void spinodal(int n_x, int n_y, double ave_comp, double *comp){

int i1, i2;

for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_x; ++i2){
comp[i2+n_y*i1] = ave_comp;
}}
	
}

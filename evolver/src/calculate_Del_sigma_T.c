
void calculate_Del_sigma_T(double DeltaC[3][3][3][3], 
double **epsilon_T, double **Del_sigma_T){

int i1,i2,i3,i4;

/* Initialise */

Del_sigma_T[1][1] = 0.0;
Del_sigma_T[1][2] = 0.0;
Del_sigma_T[2][1] = 0.0;
Del_sigma_T[2][2] = 0.0;

/* Get the eigenstress - Corresponding to DeltaC */

for(i1=1; i1<3; ++i1){
for(i2=1; i2<3; ++i2){
for(i3=1; i3<3; ++i3){
for(i4=1; i4<3; ++i4){
Del_sigma_T[i1][i2] = Del_sigma_T[i1][i2] 
	+ DeltaC[i1][i2][i3][i4]*epsilon_T[i3][i4];
}}}}

}

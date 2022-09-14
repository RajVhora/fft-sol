/** My header file for initialisation function declarations **/

extern void generate_initial_composition_profile(int n_x, int n_y, 
double *comp);

extern void add_noise(int n_x, int n_y, double *comp, double
noise_str, int SNA, int NOISEADDER);

extern void random_ppt(int n_x, int n_y, int R, int N, 
double comp_precipitate, double ave_comp, double *comp);

extern void circular(int n_x, int n_y, int R, 
double comp_precipitate, double comp_matrix, double *comp);

extern void spinodal(int n_x, int n_y, double ave_comp, double *comp);

extern void elliptic(int n_x, int n_y, double R1, double R2, 
double comp_precipitate, double comp_matrix, double *comp);

extern void thinfilm(int n_x, int n_y, int h, double comp_precipitate,
double comp_matrix, double *comp);

extern void twocirc(int n_x, int n_y, int R, 
double comp_precipitate, double comp_matrix, double *comp);


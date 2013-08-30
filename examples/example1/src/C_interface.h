#define FALSE 0
#define TRUE 1

extern void single_corrector_step_f_(double *x,double *dx, double *out_vars, int *conv);
extern void calculate_residual_f_(double *x, double *res, int *conv);

/* Pointless repetition: remove later when proper libraries invoked 
void outputPoint(double *y);
void vect_copy(double* x, int size, const double* x1);
void vect_sum(double* x, double* x1, int size, double* sum);
void vect_prod(double* x, double a, int size, double* prod); */
double norm(double* x, int size);
/* void calculateTangent(double* y, double* y_old, int size, double* Tangent);
void makePrediction(double *y_0, double *y, double *Tangent, double dlambda);
*/

double calc_mean(double *vector, int size);
void calc_mean_class(double *vector, int *class_vector, int *nclass, int size,double *store);
void calc_mean_class_I(double *vector, int *class_vector, int *nclass, int start, int end, int *ndind, double *store);
  
double calc_sd(double *vector, double mean, int size);
void calc_sd_class(double *vector, int *class_vector, int *nclass, int nsample, double *store, double *mean);
void calc_sd_class_I(double *vector, int *class_vector, int *nclass, int start, int end,int *ndind, double *store, double *mean);

double weighted_average(double *vector, int *weights, int size);

int which_max(int *vector, int size);
  
void swap_integers(int *one, int *two);

void print_array_int(int *vector, int size);
void print_array_double(double *vector, int size);
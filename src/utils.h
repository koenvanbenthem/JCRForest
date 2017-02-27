double calc_mean(double *vector, int size);
void calc_mean_class(double *vector, int *class_vector, int *nclass, int size,double *store);
void calc_mean_class_I(double *vector, int *class_vector, int *nclass, int start, int end, int *ndind, double *store);
  
double calc_sd(double *vector, double mean, int size);
void calc_sd_class(double *vector, int *class_vector, int *nclass, int nsample, double *store, double *mean);
void calc_sd_class_I(double *vector, int *class_vector, int *nclass, int start, int end,int *ndind, double *store, double *mean);

int which_max(int *vector, int size);
  
void swap_integers(int *one, int *two);
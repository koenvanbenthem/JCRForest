#include <math.h>

void swap_integers(int *one, int *two){
  int temp = *one;
  *one = *two;
  *two = temp;
}

int which_max(int *vector, int size){
  int max = vector[0];
  int best = 0;
  for(int i = 1; i < size; i++){
    if(vector[i] > max){
      best = i;
      max = vector[i];
    }
  }
  
  return best;
}

double calc_mean(double *vector, int size){
  
  double result = 0.0;
  
  for(int i = 0; i < size; i++){
    result += vector[i] / ((double) size);
  }
  
  return result;
}

void calc_mean_class(double *vector, int *class_vector, int *nclass, int nsample,double *store){
  int size[*nclass];
  for(int i=0; i<*nclass; i++){
    store[i] = 0;
    size[i] = 0;
  }
  
  for(int i=0; i<nsample; i++){
    store[class_vector[i]-1] += vector[i];
    size[class_vector[i]-1] += 1;
  }
  
  for(int i=0; i<*nclass; i++){
    store[i] = store[i] / ((double) size[i]);
  }
}

void calc_mean_class_I(double *vector, int *class_vector, int *nclass, int start, int end, int *ndind, double *store){
  int size[*nclass];
  for(int i=0; i<*nclass; i++){
    store[i] = 0;
    size[i] = 0;
  }
  
  for(int i=start; i<end; i++){
    store[class_vector[ndind[i]]-1] += vector[ndind[i]];
    size[class_vector[ndind[i]]-1] += 1;
  }
  
  for(int i=0; i < *nclass; i++){
    store[i] = store[i]/((double) size[i]);
  }
  
}

double calc_sd(double *vector, double mean, int size){
  
  double result = 0.0;
  for(int i = 0; i < size; i++){
    result += pow(vector[i]-mean,2)/((double) size - 1.0);
  }
  
  return sqrt(result);
}

void calc_sd_class(double *vector, int *class_vector, int *nclass, int nsample, double *store, double *mean){

  int size[*nclass];
  for(int i=0; i < *nclass; i++){
    size[i] = 0;
    store[i] = 0;
  }
  
  for(int i = 0; i < nsample; i++){
    store[class_vector[i]-1] += pow(vector[i]-mean[class_vector[i]-1],2);
    size[class_vector[i]-1] += 1;
  }
  
  for(int i=0; i < *nclass; i++){
    if(size[i] > 1){
      store[i] = sqrt(store[i] / ((double)(size[i]-1)));
    } else {
      store[i] = 0;
    }
  }
}

void calc_sd_class_I(double *vector, int *class_vector, int *nclass, int start, int end,int *ndind, double *store, double *mean){
  
  int size[*nclass];
  for(int i=0; i < *nclass; i++){
    size[i] = 0;
    store[i] = 0;
  }
  
  for(int i = start; i < end; i++){
    store[class_vector[ndind[i]]-1] += pow(vector[ndind[i]]-mean[class_vector[ndind[i]]-1],2);
    size[class_vector[ndind[i]]-1] += 1;
  }
  
  for(int i=0; i < *nclass; i++){
    if(size[i] > 1){
      store[i] = sqrt(store[i] / ((double)(size[i]-1)));
    } else {
      store[i] = 0;
    }
  }
}
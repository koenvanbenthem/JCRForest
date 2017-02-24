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

void calc_mean_class(double *vector, int *class_vector, int *nclass, int size,double* store){
  for(int i=0; i<*nclass; i++) store[i] = 0;
  
  for(int i=0; i<size; i++){
    store[class_vector[i]] += vector[i] / ((double) size);
  }
}
/*
void calc_mean_class_I(double *vector,){
  
}*/

double calc_sd(double *vector, double mean, int size){
  
  double result = 0.0;
  for(int i = 0; i < size; i++){
    result += pow(vector[i]-mean,2)/((double) size - 1.0);
  }
  
  return sqrt(result);
}
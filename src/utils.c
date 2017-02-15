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


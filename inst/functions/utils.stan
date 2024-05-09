
  // Create vector with elements 1, ..., N
  vector seq_len(data int N) {
    vector[N] v = rep_vector(1.0, N);
    for(n in 2:N) v[n] = n;
    return(v);
  }

    // Create integer array with elements 1, ..., N
  array[] int seq_len_int(data int N) {
    array[N] int v = rep_array(1, N);
    for(n in 2:N) v[n] = n;
    return(v);
  }

  // Same as rep(x, times=N) in R
  vector rep_vector_times(vector x, data int N) {
    return to_vector(rep_matrix(x, N));
  }

  // Repeat a matrix horizontally
  matrix rep_matrix_times(data matrix x, data int times) {
    int N = rows(x);
    int M = cols(x);
    matrix[N, M*times] out;
    for(j in 1:times) {
      out[:, (1+(j-1)*M):(j*M)] = x;
    }
    return(out);
  }

    // Find indices of val
  array[] int which(array[] int vals, int val){
    int n = size(vals);
    array[n] int out;
    int J = 0;
    for(i in 1:n){
      if(i==val){
        out[n] = i;
        J = J + 1;
      }
    }
    return(out[1:J]);
  }

  // Find first occurence of value in array
  int first(array[] int arr, int val){
    int J = size(arr);
    for(j in 1:J){
      if(arr[j]==val){
        return(j);
      }
    }
    return(0);
  }

  // Find last occurence of value in array
  int last(array[] int arr, int val){
    int J = size(arr);
    array[J] int r_arr = reverse(arr);
    for(j in 1:J){
      if(r_arr[j]==val){
        return(J-j+1);
      }
    }
    return(0);
  }

  // Create an index array where first column is start indices and second has
  // the end indices
  // Assumes that arr is sorted
  array[,] int inds_array(array[] int arr, int G){
    array[G, 2] int out;
    for(g in 1:G){
      out[g, 1] = first(arr, g);
      out[g, 2] = last(arr, g);
    }
    print("Created index array:");
    return(out);
  }

  // Transform parameter to log natural scale (vectorized)
  // - all vector arguments should have the same length
  vector to_log_natural_scale(vector log_z, vector log_mu, vector log_sigma) {
      return(log_sum_exp(log_z + log_sigma, log_mu));
  }



  // Basis function matrix (EQ kernel)
  matrix bf_eq(vector x, data matrix mat_B, data real L) {
    int N = num_elements(x);
    int B = cols(mat_B);
    matrix[N, B] mat_X = rep_matrix(x+L, B);
    matrix[N, B] PHI = 1.0/sqrt(L)*sin(0.5*pi()/L * mat_X .* mat_B);
    return(PHI);
  }

  // Compute spectral density of EQ kernel
  vector log_spd_eq(real alpha, real ell, vector omega){
    return 2*log(alpha)+log(ell)+0.5*log(2*pi())-0.5*ell^2*omega .* omega;
  }

  // Compute the multipliers s_b
  vector bf_eq_multips(real alpha, real ell, data vector seq_B,
      data real L)
  {
    return exp(0.5*log_spd_eq(alpha, ell, 0.5*pi()*seq_B/L));
  }


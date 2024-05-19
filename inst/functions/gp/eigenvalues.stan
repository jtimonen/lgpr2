  // Square roots of Eigenvalues
  vector eigval_shared(real alpha, real ell, data vector seq_B, data real L){
    return(bf_eq_multips(alpha, ell, seq_B, L));
  }

  // Square roots of Eigenvalues
  vector eigval_grouped(real alpha, real ell, data vector seq_B, data real L,
      data int G){
    int B = num_elements(seq_B);
    vector[B] s = bf_eq_multips(alpha, ell, seq_B, L);
    vector[G-1] d = bf_zs_multips(G);
    return(bf_eqXzs_multips(s, d));
  }

  // Compute spectral density of EQ kernel
  vector log_spd_eq(real alpha, real ell, vector omega){
    return 2*log(alpha)+log(ell)+0.5*log(2*pi())-0.5*ell^2*omega .* omega;
  }

  // Compute the multipliers s_b
  vector bf_eq_multips(real alpha, real ell, data vector seq_B, data real L){
    return exp(0.5*log_spd_eq(alpha, ell, 0.5*pi()*seq_B/L));
  }

  // Compute the multipliers d_g
  vector bf_zs_multips(int G){
    int Gm1 = G - 1;
    real d = sqrt(1.0+1.0/Gm1);
    return(rep_vector(d, Gm1));
  }

  // Compute the multipliers d_h * s_b
  // - s = sqrts. of eigenvalues of the basis functions (size num_bf)
  // - d = sqrts. of eigenvalues of the basis functions (size num_groups)
  vector bf_eqXzs_multips(vector s, vector d) {
    int B = num_elements(s);
    int G = num_elements(d);
    vector[B*G] delt = rep_vector(0, B*G);
    for(g in 1:G){
      delt[((g-1)*B):(g*B)] = s * d[g];
    }
    return(delt);
  }

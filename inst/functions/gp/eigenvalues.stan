  // Square roots of Eigenvalues
  vector eigval_shared(real alpha, real ell, data vector seq_B, data real L){
    return(bf_eq_multips(alpha, ell, seq_B, L));
  }

  // Square roots of Eigenvalues
  vector eigval_grouped(real alpha, real ell, data vector seq_B, data real L,
      data int G){
    int B = num_elements(seq_B);
    int Gm1 = G - 1;
    vector[B] s = bf_eq_multips(alpha, ell, seq_B, L);
    vector[B*Gm1] DELTA;
    for(g in 1:Gm1){
      DELTA[((g-1)*B+1):(g*B)] = sqrt(1.0+1.0/Gm1) * s;
    }
    return(DELTA);
  }

  // Compute spectral density of EQ kernel
  vector log_spd_eq(real alpha, real ell, vector omega){
    return 2*log(alpha)+log(ell)+0.5*log(2*pi())-0.5*ell^2*omega .* omega;
  }

  // Compute the multipliers s_b
  vector bf_eq_multips(real alpha, real ell, data vector seq_B, data real L){
    return exp(0.5*log_spd_eq(alpha, ell, 0.5*pi()*seq_B/L));
  }

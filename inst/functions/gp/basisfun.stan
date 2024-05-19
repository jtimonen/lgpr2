
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


  // Basis function matrix (zerosum kernel)
  matrix bf_zs(array[] int z, int G) {
    int N = size(z);
    int Gm1 = G - 1;
    matrix[N, G-1] VARPHI;
    matrix[G, G] H = ones_cbind_helmert_norm(G);
    print("H = ", H);
    for(g in 2:G){
      VARPHI[:, g-1] = H[z, g];
    }
    return(VARPHI);
  }

  // Compute the multipliers d_g
  vector bf_zs_multips(int G){
    int Gm1 = G - 1;
    real d = sqrt(1.0+1.0/Gm1);
    return(rep_vector(d, Gm1));
  }

  // Helmert contrast matrix (with column of ones added as first column)
  matrix ones_cbind_helmert(int G){
    matrix[G, G] H;
    for(r in 1:G){
      for(c in 1:G){
        if(c == 1){
          H[r,c] = 1;
        } else {
          if(r < c){
            H[r,c] = -1;
          } else if (r==c) {
            H[r,c] = c-1;
          } else {
            H[r,c] = 0;
          }
        }
      }
    }
    return(H);
  }


  // Helmert contrast matrix with normalized rows (with col of 1s as 1st col)
  matrix ones_cbind_helmert_norm(int G){
    matrix[G, G] H = ones_cbind_helmert(G);
    for(g in 2:G){
      H[:,g] = H[:,g]/norm2(H[:,g]);
    }
    return(H);
  }

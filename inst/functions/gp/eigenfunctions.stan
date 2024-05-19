  // Eigen functions for shared term
  matrix eigfun_shared(vector x, data matrix mat_B, data real L){
    return(bf_eq(x, mat_B, L));
  }

  // Eigen functions for grouped term
  matrix eigfun_grouped(vector x, data matrix mat_B, data real L,
      data array[] int z, data int G){
    int N = num_elements(x);
    int B = cols(mat_B);
    matrix[N, B] PHI = bf_eq(x, mat_B, L);
    matrix[N, G-1] VARPHI = bf_zs(z, G);
    return(bf_eqXzs(PHI, VARPHI));
  }

  // Basis function matrix (EQ kernel)
  matrix bf_eq(vector x, data matrix mat_B, data real L) {
    int N = num_elements(x);
    int B = cols(mat_B);
    matrix[N, B] mat_X = rep_matrix(x+L, B);
    matrix[N, B] PHI = 1.0/sqrt(L)*sin(0.5*pi()/L * mat_X .* mat_B);
    return(PHI);
  }

  // Basis function matrix (zerosum kernel)
  matrix bf_zs(data array[] int z, data int G) {
    int N = size(z);
    int Gm1 = G - 1;
    matrix[N, G-1] VARPHI;
    matrix[G, G-1] H = helmert_norm(G);
    for(g in 2:G){
      VARPHI[:, g-1] = H[z, g-1];
    }
    return(VARPHI);
  }

  // Compute product basis functions (interaction)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - VARPHI = the evaluated cat basis functions (size num_obs x (G-1))
  matrix bf_eqXzs(matrix PHI, matrix VARPHI) {
    int N = rows(PHI);
    int B = cols(PHI);
    int Gm1 = cols(VARPHI);
    matrix[N, B*Gm1] PSI;
    int idx;
    for(g in 1:Gm1){
      for(b in 1:B){
        idx = (g-1)*B + b;
        PSI[:, idx] = PHI[:,b] .* VARPHI[:, g];
      }
    }
    return(PSI);
  }



  // Helmert contrast matrix
  matrix helmert(int G){
    matrix[G, G-1] H;
    for(r in 1:G){
      for(c in 2:G){
        if(r < c){
          H[r,c-1] = -1;
        } else if (r==c) {
          H[r,c-1] = c-1;
        } else {
          H[r,c-1] = 0;
        }
      }
    }
    return(H);
  }

  // Helmert contrast matrix with normalized columns
  matrix helmert_norm(int G){
    matrix[G, G-1] H = helmert(G);
    for(g in 2:G){
      H[:,g-1] = H[:,g-1]/norm2(H[:,g-1]);
    }
    return(H);
  }

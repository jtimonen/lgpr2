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
    matrix[G, G] H = ones_cbind_helmert_norm(G);
    for(g in 2:G){
      VARPHI[:, g-1] = H[z, g];
    }
    return(VARPHI);
  }

  // Compute product basis functions (interaction)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - VARPHI = the evaluated cat basis functions (size num_obs x num_groups)
  matrix bf_eqXzs(matrix PHI, matrix VARPHI) {
    int N = rows(PHI);
    int B = cols(PHI);
    int Gm1 = cols(VARPHI)-1;
    matrix[N, B*Gm1] PSI;
    int idx = 0;
    for(g in 1:Gm1){
      for(b in 1:B){
        idx = (g-1)*B + b;
        PSI[:, idx] = PHI[:,b] .* VARPHI[:, g];
      }
    }
    return(PSI);
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

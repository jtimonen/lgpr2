  // Helmert contrast matrix
  matrix helmert(int G){
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

  // Compute a group-specific GP effect (interaction)
  // - xi = the basis function coefficient parameters (size num_groups-1 x num_bf)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - s = sqrts. of eigenvalues of the basis functions (size num_bf)
  // - group = group assignment of each data point (size num_obs)
  vector compute_f_group(matrix xi, matrix PHI, vector s, array[] int group) {
    int Gm1 = rows(xi);
    int B = num_elements(s);
    int N = rows(PHI);
    vector[Gm1] d;
    vector[N] f = rep_vector(0, N);
    matrix[N, Gm1] VARPHI;
    matrix[Gm1+1, Gm1+1] H = helmert(Gm1+1);
    for(g in 1:Gm1){
      d[g] = sqrt(1+1/Gm1);
      VARPHI[:, g] = H[group, g+1];
    }
    for(g in 1:Gm1){
      for(b in 1:B){
        f = f + xi[g,b] * s[b] * d[g] * PHI[:,b] .* VARPHI[:, g];
      }
    }
    return(f);
  }

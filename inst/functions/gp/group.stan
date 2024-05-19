  // Compute a group-specific GP effect (interaction)
  // - xi = the basis function coefficient parameters (size num_groups-1 x num_bf)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - s = sqrts. of eigenvalues of the basis functions (size num_bf)
  // - VARPHI = the evaluated cat basis functions (size num_obs x num_groups)
  // - d = sqrts. of eigenvalues of the basis functions (size num_groups)
  vector compute_f_group(matrix xi, matrix PHI, vector s,
      matrix VARPHI, vector d) {
    int Gm1 = rows(xi);
    int B = num_elements(s);
    int N = rows(PHI);
    vector[N] f = rep_vector(0, N);

    for(g in 1:Gm1){
      for(b in 1:B){
        f = f + xi[g,b] * s[b] * d[g] * PHI[:,b] .* VARPHI[:, g];
      }
    }
    return(f);
  }




  // Compute a group-specific GP effect (interaction)
  // - xi = the basis function coefficient parameters (size num_groups x num_bf)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - s = sqrts. of eigenvalues of the basis functions (size num_bf)
  // - group = group assignment of each data point (size num_obs)
  vector compute_f_group(matrix xi, matrix PHI, vector s, array[] int group) {
    int num_obs = rows(PHI);
    vector[num_obs] f;
    for(n in 1:num_obs){
      f[n] = sum(to_vector(PHI[n,:]) .* s .* to_vector(xi[group[n], :]));
    }
    return(f);
  }

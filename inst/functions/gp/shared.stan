
  // Compute a shared GP effect
  // - xi = the basis function coefficient parameters (size num_bf)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - s = sqrts. of eigenvalues of the basis functions (size num_bf)
  vector compute_f_shared(vector xi, matrix PHI, vector s) {
    return(PHI * (s .* xi));
  }

functions {
  // The usual integrand:
  // D_C = c/H0 int(cosmo_integrand(z), z = 0..z)
  real[] cosmo_integrand(real z, real[] Ez, real[] cosmo_params, real[] data_reals, int[] data_ints) {
    real Om;
    real Ol;
    real zp1;
    real H0;

    real result[1];
    
    Om = cosmo_params[1];
    Ol = cosmo_params[2];

    zp1 = z + 1.0;

    result[1] = 1.0/sqrt(Om*zp1^3 + Ol); // Enforce flat universe, so ignore Omegak
    return result;    
  }

  real[] dls(int Nobs, real[] zs, real H0, real[] Omegas, real[] data_reals, int[] data_ints) {
    real dm_matrix[Nobs, 1]; // The ODE solvers output a matrix of
			     // states at each time; in our example,
			     // the state in one-dimensional, so the
			     // matrix has only a single column
    real dl[Nobs];

    real dl0[1];

    dl0[1] = 0.0;
    
    dm_matrix = integrate_ode_rk45(cosmo_integrand, dl0, 0.0, zs, Omegas, data_reals, data_ints);

    for (i in 1:Nobs) {
      dl[i] = (1.0 + zs[i]) * 2.99792e5 / H0 * dm_matrix[i,1];
    }

    return dl;
  }
}

data {
  int Nobs;
  real dls_obs[Nobs]; // In Mpc
  real sigma_dls_obs[Nobs]; // In Mpc
  real zobs[Nobs]; // Assume redshifts are known perfectly.
}

transformed data {
  // We don't need any data for the ODE integrator.
  real data_reals[0];
  int data_ints[0];
}

parameters {
  real<lower=0> H0; // In km/s/Mpc!
  simplex[2] Omegas; // Must sum to one!
}

transformed parameters {
  real dls_pred[Nobs];

  dls_pred = dls(Nobs, zobs, H0, to_array_1d(Omegas), data_reals, data_ints);
}

model {
  dls_obs ~ normal(dls_pred, sigma_dls_obs);

  H0 ~ normal(70, 10);
  // Flat prior on Omegas
}

functions {
  // The usual integrand:
  // D_C = c/H0 int(cosmo_integrand(z), z = 0..z)
  real[] cosmo_integrand(real z, real[] Ez, real[] cosmo_params, real[] data_reals, int[] data_ints) {
    real Om;
    real Ol;
    real Ok;
    real zp1;

    real result[1];

    Om = cosmo_params[1];
    Ol = cosmo_params[2];
    Ok = 1.0 - Om - Ol;

    zp1 = z + 1.0;

    result[1] = 1.0/sqrt(Om*zp1^3 + Ol + Ok*zp1^2);
    return result;
  }

  real[] dls(real[] zs, real H0, real[] Omegas, real[] data_reals, int[] data_ints) {
    int Nobs = size(zs);
    real dm_matrix[Nobs, 1]; // The ODE solvers output a matrix of
			     // states at each time; in our example,
			     // the state in one-dimensional, so the
			     // matrix has only a single column
    real dl[Nobs];

    real dm0[1];

    dm0[1] = 0.0;

    dm_matrix = integrate_ode_rk45(cosmo_integrand, dm0, 0.0, zs, Omegas, data_reals, data_ints);

    for (i in 1:Nobs) {
      real Ok = 1.0 - sum(Omegas);
      real dM;

      if (Ok > 0) {
        dM = 1.0 / sqrt(Ok) * sinh(sqrt(Ok)*dm_matrix[i,1]);
      } else if (Ok < 0) {
        dM = 1.0 / sqrt(-Ok) * sin(sqrt(-Ok)*dm_matrix[i,1]);
      } else {
        dM = dm_matrix[i,1];
      }
      dl[i] = (1.0 + zs[i]) * 2.99792e5 / H0 * dM;
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
  real<lower=0, upper=1> Omegas[2]; // Omega_M, Omega_L
}

transformed parameters {
  real Ok = 1 - sum(Omegas);
  real dls_pred[Nobs];

  dls_pred = dls(zobs, H0, Omegas, data_reals, data_ints);
}

model {
  dls_obs ~ normal(dls_pred, sigma_dls_obs);

  H0 ~ normal(70, 10);
  // Flat prior on Omegas
}

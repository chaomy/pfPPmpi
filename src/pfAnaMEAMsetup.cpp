/*
 * @Author: chaomy
 * @Date:   2018-01-30 13:42:16
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-04 22:22:52
 */

#include "pfHome.h"

void pfHome::meam_setup_global(const arma::mat& vv) {
  int cn = 0;
  for (int i = 0; i < nelt; i++) {
    alpha_meam[i][i] = vv[cn++];
    beta0_meam[i] = vv[cn++];
    beta1_meam[i] = vv[cn++];
    beta2_meam[i] = vv[cn++];
    beta3_meam[i] = vv[cn++];
    Ec_meam[i][i] = vv[cn++];
    A_meam[i] = vv[cn++];
    t0_meam[i] = vv[cn++];
    t1_meam[i] = vv[cn++];
    t2_meam[i] = vv[cn++];
    t3_meam[i] = vv[cn++];
    rc_meam = vv[cn++];
  }
}

void pfHome::meam_setup_global(const vector<double>& vv) {
  int cn = 0;
  for (int i = 0; i < nelt; i++) {
    alpha_meam[i][i] = vv[cn++];
    beta0_meam[i] = vv[cn++];
    beta1_meam[i] = vv[cn++];
    beta2_meam[i] = vv[cn++];
    beta3_meam[i] = vv[cn++];
    Ec_meam[i][i] = vv[cn++];
    A_meam[i] = vv[cn++];
    t0_meam[i] = vv[cn++];
    t1_meam[i] = vv[cn++];
    t2_meam[i] = vv[cn++];
    t3_meam[i] = vv[cn++];
    rc_meam = vv[cn++];
  }
}

void pfHome::meam_setup_globalfixed() {  //  those are fixed
  lattp = vector<lattice_t>({BCC});
  for (int i = 0; i < nelt; i++) {
    cout << "rho0 = " << (rho0_meam[i] = rozero[i]) << endl;
    cout << "ibar = " << (ibar_meam[i] = ibar[i]) << endl;
    cout << "lat = " << (lattce_meam[i][i] = lattp[i]) << endl;
    cout << "z = " << (Z_meam[i] = cnn1[i]) << endl;
    cout << "ielt = " << (ielt_meam[i] = 1) << endl;

    if (lattce_meam[i][i] == FCC)
      re_meam[i][i] = alat[i] / sqrt(2.0);
    else if (lattce_meam[i][i] == BCC)
      re_meam[i][i] = alat[i] * sqrt(3.0) / 2.0;
    else if (lattce_meam[i][i] == HCP)
      re_meam[i][i] = alat[i];
    else if (lattce_meam[i][i] == DAM)
      re_meam[i][i] = alat[i];
    else if (lattce_meam[i][i] == DIA)
      re_meam[i][i] = alat[i] * sqrt(3.0) / 4.0;

    cout << "alat = " << alat[i] << endl;
    cout << "re = " << re_meam[i][i] << endl;
  }

  // setup boundaries alpha  b0  b1   b2   b3   Ec   A t0 t1 t2 t3 rc_meam
  lob =
      vector<double>({3, -5.0, -5.0, -5.0, -5.0, 1.0, 0.0, 0, -5, -5, -5, 4.7});
  hib = vector<double>({6, 5.0, 5.0, 5.0, 5.0, 10, 2.0, 5, 5., 5., 5., 6.25});
  for (int i = 0; i < lob.size(); i++) deb.push_back(hib[i] - lob[i]);
}

void pfHome::meam_setup_global() {  // for benchmark, and reference
  elems = vector<string>({"Mg"});
  lattp = vector<lattice_t>({HCP});
  cnn1 = vector<int>({12});
  meamparms["atwt"] = vector<double>({24.320});
  meamparms["alpha"] = vector<double>({5.5663414489});
  meamparms["b0"] = vector<double>({2.3});
  meamparms["b1"] = vector<double>({1.0});
  meamparms["b2"] = vector<double>({3.0});
  meamparms["b3"] = vector<double>({1.0});
  meamparms["alat"] = vector<double>({3.2});
  meamparms["esub"] = vector<double>({1.55});
  meamparms["asub"] = vector<double>({0.52});
  meamparms["t0"] = vector<double>({1.0});
  meamparms["t1"] = vector<double>({9.0});
  meamparms["t2"] = vector<double>({-2.0});
  meamparms["t3"] = vector<double>({-9.5});
  meamparms["rozero"] = vector<double>({1.0});
  ielement = vector<int>({1});
  ibar = vector<int>({3});

  for (int i = 0; i < nelt; i++) {
    lattce_meam[i][i] = lattp[i];
    Z_meam[i] = cnn1[i];
    ielt_meam[i] = ielement[i];
    alpha_meam[i][i] = meamparms["alpha"][i];
    beta0_meam[i] = meamparms["b0"][i];
    beta1_meam[i] = meamparms["b1"][i];
    beta2_meam[i] = meamparms["b2"][i];
    beta3_meam[i] = meamparms["b3"][i];
    Ec_meam[i][i] = meamparms["esub"][i];
    A_meam[i] = meamparms["asub"][i];
    t0_meam[i] = meamparms["t0"][i];
    t1_meam[i] = meamparms["t1"][i];
    t2_meam[i] = meamparms["t2"][i];
    t3_meam[i] = meamparms["t3"][i];
    rho0_meam[i] = meamparms["rozero"][i];
    ibar_meam[i] = ibar[i];
  }
}

void pfHome::meam_setup_done() {
  cutforce = rc_meam;
  cutforcesq = cutforce * cutforce;

  // augment t1 term
  for (int i = 0; i < nelt; i++) t1_meam[i] += augt1 * 3.0 / 5.0 * t3_meam[i];
  alloyparams();  // indices and factors for Voight notation

  // indices and factors for Voight notation
  int nv2 = 0, nv3 = 0, m, n, p;
  for (m = 0; m < 3; m++) {
    for (n = m; n < 3; n++) {
      vind2D[m][n] = nv2;
      vind2D[n][m] = nv2;
      nv2++;
      for (p = n; p < 3; p++) {
        vind3D[m][n][p] = nv3;
        vind3D[m][p][n] = nv3;
        vind3D[n][m][p] = nv3;
        vind3D[n][p][m] = nv3;
        vind3D[p][m][n] = nv3;
        vind3D[p][n][m] = nv3;
        nv3++;
      }
    }
  }

  v2D[0] = 1, v2D[1] = 2, v2D[2] = 2;
  v2D[3] = 1, v2D[4] = 2, v2D[5] = 1;

  v3D[0] = 1, v3D[1] = 3, v3D[2] = 3;
  v3D[3] = 3, v3D[4] = 6, v3D[5] = 3;
  v3D[6] = 1, v3D[7] = 3, v3D[8] = 3;
  v3D[9] = 1;

  nv2 = 0;
  for (m = 0; m < nelt; m++) {
    for (n = m; n < nelt; n++) {
      eltind[m][n] = nv2;
      eltind[n][m] = nv2;
      nv2++;
    }
  }

  // compute background densities for reference structure
  compute_reference_density();

  // compute pair potential and setup array for interpolation
  nr = 1000;
  dr = 1.1 * rc_meam / nr;
  compute_pair_meam();
}

void pfHome::alloyparams() {
  int i, j, k;
  double eb;

  // Loop over pairs
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nelt; j++) {
      // Treat off-diagonal pairs If i>j, set all equal to i<j case
      // (which has aready been set,  here or in the input file)
      if (i > j) {
        re_meam[i][j] = re_meam[j][i];
        Ec_meam[i][j] = Ec_meam[j][i];
        alpha_meam[i][j] = alpha_meam[j][i];
        lattce_meam[i][j] = lattce_meam[j][i];
        nn2_meam[i][j] = nn2_meam[j][i];
        // If i<j and term is unset, use default values (e.g. mean of i-i and
        // j-j)
      } else if (j > i) {
        if (iszero(Ec_meam[i][j])) {
          if (lattce_meam[i][j] == L12)
            Ec_meam[i][j] =
                (3 * Ec_meam[i][i] + Ec_meam[j][j]) / 4.0 - delta_meam[i][j];
          else if (lattce_meam[i][j] == C11) {
            if (lattce_meam[i][i] == DIA)
              Ec_meam[i][j] =
                  (2 * Ec_meam[i][i] + Ec_meam[j][j]) / 3.0 - delta_meam[i][j];
            else
              Ec_meam[i][j] =
                  (Ec_meam[i][i] + 2 * Ec_meam[j][j]) / 3.0 - delta_meam[i][j];
          } else
            Ec_meam[i][j] =
                (Ec_meam[i][i] + Ec_meam[j][j]) / 2.0 - delta_meam[i][j];
        }
        if (iszero(alpha_meam[i][j]))
          alpha_meam[i][j] = (alpha_meam[i][i] + alpha_meam[j][j]) / 2.0;
        if (iszero(re_meam[i][j]))
          re_meam[i][j] = (re_meam[i][i] + re_meam[j][j]) / 2.0;
      }
    }
  }

  // Cmin[i][k][j] is symmetric in i-j, but not k.  For all triplets
  // where i>j, set equal to the i<j element.  Likewise for Cmax.
  for (i = 1; i < nelt; i++) {
    for (j = 0; j < i; j++) {
      for (k = 0; k < nelt; k++) {
        Cmin_meam[i][j][k] = Cmin_meam[j][i][k];
        Cmax_meam[i][j][k] = Cmax_meam[j][i][k];
      }
    }
  }

  // ebound gives the squared distance such that, for rik2 or rjk2>ebound,
  // atom k definitely lies outside the screening function ellipse (so
  // there is no need to calculate its effects).  Here, compute it for all
  // triplets [i][j][k] so that ebound[i][j] is the maximized over k
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nelt; j++) {
      for (k = 0; k < nelt; k++) {
        eb = (Cmax_meam[i][j][k] * Cmax_meam[i][j][k]) /
             (4.0 * (Cmax_meam[i][j][k] - 1.0));
        ebound_meam[i][j] = std::max(ebound_meam[i][j], eb);
      }
    }
  }
}

//----------------------------------------------------------------------c
// Compute MEAM pair potential for distance r, element types a and b
//
double pfHome::phi_meam(double r, int a, int b) {
  /*unused:double a1,a2,a12;*/
  double t11av, t21av, t31av, t12av, t22av, t32av;
  double G1, G2, s1[3], s2[3], rho0_1, rho0_2;
  double Gam1, Gam2, Z1, Z2;
  double rhobar1, rhobar2, F1, F2;
  double rho01, rho11, rho21, rho31;
  double rho02, rho12, rho22, rho32;
  double scalfac, phiaa, phibb;
  double Eu;
  double arat, scrn /*unused:,scrn2*/;
  int Z12, errorflag;
  int n, nmax, Z1nn, Z2nn;
  lattice_t latta /*unused:,lattb*/;
  double rho_bkgd1, rho_bkgd2;

  double phi_m = 0.0;

  // Equation numbers below refer to:
  //   I. Huang et.al., Modelling simul. Mater. Sci. Eng. 3:615

  // get number of neighbors in the reference structure
  //   Nref[i][j] = # of i's neighbors of type j

  Z12 = get_Zij(lattce_meam[a][b]);
  get_densref(r, a, b, &rho01, &rho11, &rho21, &rho31, &rho02, &rho12, &rho22,
              &rho32);

  if (rho01 <= 1e-14 && rho02 <= 1e-14) return 0.0;

  // calculate average weighting factors for the reference structure
  if (lattce_meam[a][b] == C11) {
    if (ialloy == 2) {
      t11av = t1_meam[a];
      t12av = t1_meam[b];
      t21av = t2_meam[a];
      t22av = t2_meam[b];
      t31av = t3_meam[a];
      t32av = t3_meam[b];
    } else {
      scalfac = 1.0 / (rho01 + rho02);
      t11av = scalfac * (t1_meam[a] * rho01 + t1_meam[b] * rho02);
      t12av = t11av;
      t21av = scalfac * (t2_meam[a] * rho01 + t2_meam[b] * rho02);
      t22av = t21av;
      t31av = scalfac * (t3_meam[a] * rho01 + t3_meam[b] * rho02);
      t32av = t31av;
    }
  } else {
    // average weighting factors for the reference structure, eqn. I.8
    get_tavref(&t11av, &t21av, &t31av, &t12av, &t22av, &t32av, t1_meam[a],
               t2_meam[a], t3_meam[a], t1_meam[b], t2_meam[b], t3_meam[b], r, a,
               b, lattce_meam[a][b]);
  }

  // for c11b structure, calculate background electron densities
  if (lattce_meam[a][b] == C11) {
    latta = lattce_meam[a][a];
    if (latta == DIA) {
      rhobar1 = pow(((Z12 / 2) * (rho02 + rho01)), 2) +
                t11av * pow((rho12 - rho11), 2) +
                t21av / 6.0 * pow(rho22 + rho21, 2) +
                121.0 / 40.0 * t31av * pow((rho32 - rho31), 2);
      rhobar1 = sqrt(rhobar1);
      rhobar2 = pow(Z12 * rho01, 2) + 2.0 / 3.0 * t21av * pow(rho21, 2);
      rhobar2 = sqrt(rhobar2);
    } else {
      rhobar2 = pow(((Z12 / 2) * (rho01 + rho02)), 2) +
                t12av * pow((rho11 - rho12), 2) +
                t22av / 6.0 * pow(rho21 + rho22, 2) +
                121.0 / 40.0 * t32av * pow((rho31 - rho32), 2);
      rhobar2 = sqrt(rhobar2);
      rhobar1 = pow(Z12 * rho02, 2) + 2.0 / 3.0 * t22av * pow(rho22, 2);
      rhobar1 = sqrt(rhobar1);
    }
  } else {
    // for other structures, use formalism developed in Huang's paper
    //
    //     composition-dependent scaling, equation I.7
    //     If using mixing rule for t, apply to reference structure; else
    //     use precomputed values
    if (mix_ref_t == 1) {
      Z1 = Z_meam[a];
      Z2 = Z_meam[b];
      if (ibar_meam[a] <= 0)
        G1 = 1.0;
      else {
        get_shpfcn(lattce_meam[a][a], s1);
        Gam1 = (s1[0] * t11av + s1[1] * t21av + s1[2] * t31av) / (Z1 * Z1);
        G1 = G_gam(Gam1, ibar_meam[a], errorflag);
      }
      if (ibar_meam[b] <= 0)
        G2 = 1.0;
      else {
        get_shpfcn(lattce_meam[b][b], s2);
        Gam2 = (s2[0] * t12av + s2[1] * t22av + s2[2] * t32av) / (Z2 * Z2);
        G2 = G_gam(Gam2, ibar_meam[b], errorflag);
      }
      rho0_1 = rho0_meam[a] * Z1 * G1;
      rho0_2 = rho0_meam[b] * Z2 * G2;
    }
    Gam1 = (t11av * rho11 + t21av * rho21 + t31av * rho31);
    if (rho01 < 1.0e-14)
      Gam1 = 0.0;
    else
      Gam1 = Gam1 / (rho01 * rho01);

    Gam2 = (t12av * rho12 + t22av * rho22 + t32av * rho32);
    if (rho02 < 1.0e-14)
      Gam2 = 0.0;
    else
      Gam2 = Gam2 / (rho02 * rho02);

    G1 = G_gam(Gam1, ibar_meam[a], errorflag);
    G2 = G_gam(Gam2, ibar_meam[b], errorflag);
    if (mix_ref_t == 1) {
      rho_bkgd1 = rho0_1;
      rho_bkgd2 = rho0_2;
    } else {
      if (bkgd_dyn == 1) {
        rho_bkgd1 = rho0_meam[a] * Z_meam[a];
        rho_bkgd2 = rho0_meam[b] * Z_meam[b];
      } else {
        rho_bkgd1 = rho_ref_meam[a];
        rho_bkgd2 = rho_ref_meam[b];
      }
    }
    rhobar1 = rho01 / rho_bkgd1 * G1;
    rhobar2 = rho02 / rho_bkgd2 * G2;
  }

  // compute embedding functions, eqn I.5
  if (iszero(rhobar1))
    F1 = 0.0;
  else {
    if (emb_lin_neg == 1 && rhobar1 <= 0)
      F1 = -A_meam[a] * Ec_meam[a][a] * rhobar1;
    else
      F1 = A_meam[a] * Ec_meam[a][a] * rhobar1 * log(rhobar1);
  }
  if (iszero(rhobar2))
    F2 = 0.0;
  else {
    if (emb_lin_neg == 1 && rhobar2 <= 0)
      F2 = -A_meam[b] * Ec_meam[b][b] * rhobar2;
    else
      F2 = A_meam[b] * Ec_meam[b][b] * rhobar2 * log(rhobar2);
  }
  // compute Rose function, I.16
  Eu = erose(r, re_meam[a][b], alpha_meam[a][b], Ec_meam[a][b],
             repuls_meam[a][b], attrac_meam[a][b], erose_form);

  // calculate the pair energy
  if (lattce_meam[a][b] == C11) {
    latta = lattce_meam[a][a];
    if (latta == DIA) {
      phiaa = phi_meam(r, a, a);
      phi_m = (3 * Eu - F2 - 2 * F1 - 5 * phiaa) / Z12;
    } else {
      phibb = phi_meam(r, b, b);
      phi_m = (3 * Eu - F1 - 2 * F2 - 5 * phibb) / Z12;
    }
  } else if (lattce_meam[a][b] == L12) {
    phiaa = phi_meam(r, a, a);
    //       account for second neighbor a-a potential here...
    Z1nn = get_Zij(lattce_meam[a][a]);
    Z2nn = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a], Cmax_meam[a][a][a],
                    arat, scrn);
    nmax = 10;
    if (scrn > 0.0) {
      for (n = 1; n <= nmax; n++) {
        phiaa = phiaa + pow((-Z2nn * scrn / Z1nn), n) *
                            phi_meam(r * pow(arat, n), a, a);
      }
    }
    phi_m = Eu / 3.0 - F1 / 4.0 - F2 / 12.0 - phiaa;

  } else {
    //
    // potential is computed from Rose function and embedding energy
    phi_m = (2 * Eu - F1 - F2) / Z12;
    //
  }

  // if r = 0, just return 0
  if (iszero(r)) phi_m = 0.0;

  return phi_m;
}

void pfHome::compute_reference_density() {
  int a, Z, Z2, errorflag;
  double gam, Gbar, shp[3];
  double rho0, rho0_2nn, arat, scrn;

  for (a = 0; a < nelt; a++) {  // loop over element types
    Z = (int)Z_meam[a];
    if (ibar_meam[a] <= 0)
      Gbar = 1.0;
    else {
      get_shpfcn(lattce_meam[a][a], shp);
      gam = (t1_meam[a] * shp[0] + t2_meam[a] * shp[1] + t3_meam[a] * shp[2]) /
            (Z * Z);
      Gbar = G_gam(gam, ibar_meam[a], errorflag);
    }

    //     The zeroth order density in the reference structure, with
    //     equilibrium spacing, is just the number of first neighbors times
    //     the rho0_meam coefficient...
    rho0 = rho0_meam[a] * Z;

    //     ...unless we have unscreened second neighbors, in which case we
    //     add on the contribution from those (accounting for partial
    //     screening)
    if (nn2_meam[a][a] == 1) {
      Z2 = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a], Cmax_meam[a][a][a],
                    arat, scrn);
      rho0_2nn =
          rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * (arat - 1));
      rho0 = rho0 + Z2 * rho0_2nn * scrn;
    }

    rho_ref_meam[a] = rho0 * Gbar;
  }
}

//----------------------------------------------------------------------
// Compute background density for reference structure of each element
void pfHome::compute_pair_meam() {
  double r /*ununsed:, temp*/;
  int j, a, b, nv2;
  double astar, frac, phizbl;
  int n, nmax, Z1, Z2;
  double arat, rarat, scrn, scrn2;
  double phiaa, phibb /*unused:,phitmp*/;
  double C, s111, s112, s221, S11, S22;

  phir = vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar1 =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar2 =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar3 =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar4 =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar5 =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));
  phirar6 =
      vector<vector<double>>((nelt * (nelt + 1)) / 2, vector<double>(nr, 0));

  // loop over pairs of element types
  nv2 = 0;
  for (a = 0; a < nelt; a++) {
    for (b = a; b < nelt; b++) {
      for (j = 0; j < nr; j++) {  // loop over r values and compute
        r = j * dr;
        phir[nv2][j] = phi_meam(r, a, b);

        // if using second-nearest neighbor, solve recursive problem
        // (see Lee and Baskes, PRB 62(13):8564 eqn.(21))
        if (nn2_meam[a][b] == 1) {
          Z1 = get_Zij(lattce_meam[a][b]);
          Z2 = get_Zij2(lattce_meam[a][b], Cmin_meam[a][a][b],
                        Cmax_meam[a][a][b], arat, scrn);

          //     The B1, B2,  and L12 cases with NN2 have a trick to them; we
          //     need to
          //     compute the contributions from second nearest neighbors, like
          //     a-a
          //     pairs, but need to include NN2 contributions to those pairs as
          //     well.
          if (lattce_meam[a][b] == B1 || lattce_meam[a][b] == B2 ||
              lattce_meam[a][b] == L12) {
            rarat = r * arat;
            //               phi_aa
            phiaa = phi_meam(rarat, a, a);
            Z1 = get_Zij(lattce_meam[a][a]);
            Z2 = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a],
                          Cmax_meam[a][a][a], arat, scrn);
            nmax = 10;
            if (scrn > 0.0) {
              for (n = 1; n <= nmax; n++) {
                phiaa = phiaa + pow((-Z2 * scrn / Z1), n) *
                                    phi_meam(rarat * pow(arat, n), a, a);
              }
            }

            //               phi_bb
            phibb = phi_meam(rarat, b, b);
            Z1 = get_Zij(lattce_meam[b][b]);
            Z2 = get_Zij2(lattce_meam[b][b], Cmin_meam[b][b][b],
                          Cmax_meam[b][b][b], arat, scrn);
            nmax = 10;
            if (scrn > 0.0) {
              for (n = 1; n <= nmax; n++) {
                phibb = phibb + pow((-Z2 * scrn / Z1), n) *
                                    phi_meam(rarat * pow(arat, n), b, b);
              }
            }

            if (lattce_meam[a][b] == B1 || lattce_meam[a][b] == B2) {
              //     Add contributions to the B1 or B2 potential
              Z1 = get_Zij(lattce_meam[a][b]);
              Z2 = get_Zij2(lattce_meam[a][b], Cmin_meam[a][a][b],
                            Cmax_meam[a][a][b], arat, scrn);
              phir[nv2][j] = phir[nv2][j] - Z2 * scrn / (2 * Z1) * phiaa;
              Z2 = get_Zij2(lattce_meam[a][b], Cmin_meam[b][b][a],
                            Cmax_meam[b][b][a], arat, scrn2);
              phir[nv2][j] = phir[nv2][j] - Z2 * scrn2 / (2 * Z1) * phibb;

            } else if (lattce_meam[a][b] == L12) {
              //     The L12 case has one last trick; we have to be careful to
              //     compute
              //     the correct screening between 2nd-neighbor pairs.  1-1
              //     second-neighbor pairs are screened by 2 type 1 atoms and
              //     two type
              //     2 atoms.  2-2 second-neighbor pairs are screened by 4 type
              //     1
              //     atoms.
              C = 1.0;
              get_sijk(C, a, a, a, &s111);
              get_sijk(C, a, a, b, &s112);
              get_sijk(C, b, b, a, &s221);
              S11 = s111 * s111 * s112 * s112;
              S22 = pow(s221, 4);
              phir[nv2][j] =
                  phir[nv2][j] - 0.75 * S11 * phiaa - 0.25 * S22 * phibb;
            }

          } else {
            nmax = 10;
            for (n = 1; n <= nmax; n++) {
              phir[nv2][j] =
                  phir[nv2][j] +
                  pow((-Z2 * scrn / Z1), n) * phi_meam(r * pow(arat, n), a, b);
            }
          }
        }

        // For Zbl potential:
        // if astar <= -3
        //   potential is zbl potential
        // else if -3 < astar < -1
        //   potential is linear combination with zbl potential
        // endif
        if (zbl_meam[a][b] == 1) {
          astar = alpha_meam[a][b] * (r / re_meam[a][b] - 1.0);
          if (astar <= -3.0)
            phir[nv2][j] = zbl(r, ielt_meam[a], ielt_meam[b]);
          else if (astar > -3.0 && astar < -1.0) {
            frac = fcut(1 - (astar + 1.0) / (-3.0 + 1.0));
            phizbl = zbl(r, ielt_meam[a], ielt_meam[b]);
            phir[nv2][j] = frac * phir[nv2][j] + (1 - frac) * phizbl;
          }
        }
      }
      // call interpolation
      interpolate_meam(nv2);

      nv2 = nv2 + 1;
    }  // b
  }    // a
}

//------------------------------------------------------------------------------c
// Average weighting factors for the reference structure
void pfHome::get_tavref(double* t11av, double* t21av, double* t31av,
                        double* t12av, double* t22av, double* t32av, double t11,
                        double t21, double t31, double t12, double t22,
                        double t32, double r, int a, int b, lattice_t latt) {
  double rhoa01, rhoa02, a1, a2, rho01; /*,rho02*/

  //     For ialloy = 2, no averaging is done
  if (ialloy == 2) {
    *t11av = t11;
    *t21av = t21;
    *t31av = t31;
    *t12av = t12;
    *t22av = t22;
    *t32av = t32;
  } else {
    if (latt == FCC || latt == BCC || latt == DIA || latt == HCP ||
        latt == B1 || latt == DIM || latt == B2) {
      //     all neighbors are of the opposite type
      *t11av = t12;
      *t21av = t22;
      *t31av = t32;
      *t12av = t11;
      *t22av = t21;
      *t32av = t31;
    } else {
      a1 = r / re_meam[a][a] - 1.0;
      a2 = r / re_meam[b][b] - 1.0;
      rhoa01 = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
      rhoa02 = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);
      if (latt == L12) {
        rho01 = 8 * rhoa01 + 4 * rhoa02;
        *t11av = (8 * t11 * rhoa01 + 4 * t12 * rhoa02) / rho01;
        *t12av = t11;
        *t21av = (8 * t21 * rhoa01 + 4 * t22 * rhoa02) / rho01;
        *t22av = t21;
        *t31av = (8 * t31 * rhoa01 + 4 * t32 * rhoa02) / rho01;
        *t32av = t31;
      } else {
        //      call error('Lattice not defined in get_tavref.')
      }
    }
  }
}

//------------------------------------------------------------------------------c
void pfHome::get_sijk(double C, int i, int j, int k, double* sijk) {
  double x =
      (C - Cmin_meam[i][j][k]) / (Cmax_meam[i][j][k] - Cmin_meam[i][j][k]);
  *sijk = fcut(x);
}

//------------------------------------------------------------------------------c
// Calculate density functions, assuming reference configuration
void pfHome::get_densref(double r, int a, int b, double* rho01, double* rho11,
                         double* rho21, double* rho31, double* rho02,
                         double* rho12, double* rho22, double* rho32) {
  double a1, a2;
  double s[3];
  lattice_t lat;
  int Zij2nn;
  double rhoa01nn, rhoa02nn;
  double rhoa01, rhoa11, rhoa21, rhoa31;
  double rhoa02, rhoa12, rhoa22, rhoa32;
  double arat, scrn, denom;
  double C, s111, s112, s221, S11, S22;

  a1 = r / re_meam[a][a] - 1.0;
  a2 = r / re_meam[b][b] - 1.0;

  rhoa01 = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
  rhoa11 = rho0_meam[a] * MathSpecial::fm_exp(-beta1_meam[a] * a1);
  rhoa21 = rho0_meam[a] * MathSpecial::fm_exp(-beta2_meam[a] * a1);
  rhoa31 = rho0_meam[a] * MathSpecial::fm_exp(-beta3_meam[a] * a1);
  rhoa02 = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);
  rhoa12 = rho0_meam[b] * MathSpecial::fm_exp(-beta1_meam[b] * a2);
  rhoa22 = rho0_meam[b] * MathSpecial::fm_exp(-beta2_meam[b] * a2);
  rhoa32 = rho0_meam[b] * MathSpecial::fm_exp(-beta3_meam[b] * a2);

  lat = lattce_meam[a][b];

  *rho11 = 0.0;
  *rho21 = 0.0;
  *rho31 = 0.0;
  *rho12 = 0.0;
  *rho22 = 0.0;
  *rho32 = 0.0;

  if (lat == FCC) {
    *rho01 = 12.0 * rhoa02;
    *rho02 = 12.0 * rhoa01;
  } else if (lat == BCC) {
    *rho01 = 8.0 * rhoa02;
    *rho02 = 8.0 * rhoa01;
  } else if (lat == B1) {
    *rho01 = 6.0 * rhoa02;
    *rho02 = 6.0 * rhoa01;
  } else if (lat == DIA) {
    *rho01 = 4.0 * rhoa02;
    *rho02 = 4.0 * rhoa01;
    *rho31 = 32.0 / 9.0 * rhoa32 * rhoa32;
    *rho32 = 32.0 / 9.0 * rhoa31 * rhoa31;
  } else if (lat == HCP) {
    *rho01 = 12 * rhoa02;
    *rho02 = 12 * rhoa01;
    *rho31 = 1.0 / 3.0 * rhoa32 * rhoa32;
    *rho32 = 1.0 / 3.0 * rhoa31 * rhoa31;
  } else if (lat == DAM) {
    get_shpfcn(DAM, s);
    *rho01 = rhoa02;
    *rho02 = rhoa01;
    *rho11 = s[0] * rhoa12 * rhoa12;
    *rho12 = s[0] * rhoa11 * rhoa11;
    *rho21 = s[1] * rhoa22 * rhoa22;
    *rho22 = s[1] * rhoa21 * rhoa21;
    *rho31 = s[2] * rhoa32 * rhoa32;
    *rho32 = s[2] * rhoa31 * rhoa31;
  } else if (lat == C11) {
    *rho01 = rhoa01;
    *rho02 = rhoa02;
    *rho11 = rhoa11;
    *rho12 = rhoa12;
    *rho21 = rhoa21;
    *rho22 = rhoa22;
    *rho31 = rhoa31;
    *rho32 = rhoa32;
  } else if (lat == L12) {
    *rho01 = 8 * rhoa01 + 4 * rhoa02;
    *rho02 = 12 * rhoa01;
    if (ialloy == 1) {
      *rho21 = 8. / 3. * pow(rhoa21 * t2_meam[a] - rhoa22 * t2_meam[b], 2);
      denom = 8 * rhoa01 * pow(t2_meam[a], 2) + 4 * rhoa02 * pow(t2_meam[b], 2);
      if (denom > 0.) *rho21 = *rho21 / denom * *rho01;
    } else
      *rho21 = 8. / 3. * (rhoa21 - rhoa22) * (rhoa21 - rhoa22);
  } else if (lat == B2) {
    *rho01 = 8.0 * rhoa02;
    *rho02 = 8.0 * rhoa01;
  } else {
    //        call error('Lattice not defined in get_densref.')
  }

  if (nn2_meam[a][b] == 1) {
    Zij2nn = get_Zij2(lat, Cmin_meam[a][a][b], Cmax_meam[a][a][b], arat, scrn);

    a1 = arat * r / re_meam[a][a] - 1.0;
    a2 = arat * r / re_meam[b][b] - 1.0;

    rhoa01nn = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
    rhoa02nn = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);

    if (lat == L12) {
      //     As usual, L12 thinks it's special; we need to be careful computing
      //     the screening functions
      C = 1.0;
      get_sijk(C, a, a, a, &s111);
      get_sijk(C, a, a, b, &s112);
      get_sijk(C, b, b, a, &s221);
      S11 = s111 * s111 * s112 * s112;
      S22 = pow(s221, 4);
      *rho01 = *rho01 + 6 * S11 * rhoa01nn;
      *rho02 = *rho02 + 6 * S22 * rhoa02nn;

    } else {
      //     For other cases, assume that second neighbor is of same type,
      //     first neighbor may be of different type

      *rho01 = *rho01 + Zij2nn * scrn * rhoa01nn;

      //     Assume Zij2nn and arat don't depend on order, but scrn might
      Zij2nn =
          get_Zij2(lat, Cmin_meam[b][b][a], Cmax_meam[b][b][a], arat, scrn);
      *rho02 = *rho02 + Zij2nn * scrn * rhoa02nn;
    }
  }
}

void pfHome::interpolate_meam(int ind) {
  int j;
  double drar;

  // map to coefficient space
  nrar = nr;
  drar = dr;
  rdrar = 1.0 / drar;

  // phir interp
  for (j = 0; j < nrar; j++) {
    phirar[ind][j] = phir[ind][j];
  }
  phirar1[ind][0] = phirar[ind][1] - phirar[ind][0];
  phirar1[ind][1] = 0.5 * (phirar[ind][2] - phirar[ind][0]);
  phirar1[ind][nrar - 2] =
      0.5 * (phirar[ind][nrar - 1] - phirar[ind][nrar - 3]);
  phirar1[ind][nrar - 1] = 0.0;
  for (j = 2; j < nrar - 2; j++) {
    phirar1[ind][j] = ((phirar[ind][j - 2] - phirar[ind][j + 2]) +
                       8.0 * (phirar[ind][j + 1] - phirar[ind][j - 1])) /
                      12.;
  }

  for (j = 0; j < nrar - 1; j++) {
    phirar2[ind][j] = 3.0 * (phirar[ind][j + 1] - phirar[ind][j]) -
                      2.0 * phirar1[ind][j] - phirar1[ind][j + 1];
    phirar3[ind][j] = phirar1[ind][j] + phirar1[ind][j + 1] -
                      2.0 * (phirar[ind][j + 1] - phirar[ind][j]);
  }
  phirar2[ind][nrar - 1] = 0.0;
  phirar3[ind][nrar - 1] = 0.0;

  for (j = 0; j < nrar; j++) {
    phirar4[ind][j] = phirar1[ind][j] / drar;
    phirar5[ind][j] = 2.0 * phirar2[ind][j] / drar;
    phirar6[ind][j] = 3.0 * phirar3[ind][j] / drar;
  }
}

//---------------------------------------------------------------------
// Compute Rose energy function, I.16
//
double pfHome::compute_phi(double rij, int elti, int eltj) {
  double pp;
  int ind, kk;

  ind = eltind[elti][eltj];
  pp = rij * rdrar;
  kk = (int)pp;
  kk = std::min(kk, nrar - 2);
  pp = pp - kk;
  pp = std::min(pp, 1.0);
  double result =
      ((phirar3[ind][kk] * pp + phirar2[ind][kk]) * pp + phirar1[ind][kk]) *
          pp +
      phirar[ind][kk];

  return result;
}

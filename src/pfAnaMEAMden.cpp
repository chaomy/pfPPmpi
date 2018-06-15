/*
 * @Author: chaomy
 * @Date:   2018-01-29 23:45:39
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 01:46:12
 */

#include "pfMEAMC.h"

using std::ofstream;
using std::setprecision;
using std::setw;

void pfHome::pfForce::pfMEAMC::meam_dens_setup(Config& cc) {
  int natoms = cc.natoms;

  rho = vector<double>(natoms, 0);
  rho0 = vector<double>(natoms, 0);
  rho1 = vector<double>(natoms, 0);
  rho2 = vector<double>(natoms, 0);
  rho3 = vector<double>(natoms, 0);
  frhop = vector<double>(natoms, 0);
  gamma = vector<double>(natoms, 0);
  dgamma1 = vector<double>(natoms, 0);
  dgamma2 = vector<double>(natoms, 0);
  dgamma3 = vector<double>(natoms, 0);

  arho2b = vector<double>(natoms, 0);
  arho1 = vector<vector<double>>(natoms, vector<double>(3, 0));
  arho2 = vector<vector<double>>(natoms, vector<double>(6, 0));
  arho3 = vector<vector<double>>(natoms, vector<double>(10, 0));
  arho3b = vector<vector<double>>(natoms, vector<double>(3, 0));
  t_ave = vector<vector<double>>(natoms, vector<double>(3, 0));
  tsq_ave = vector<vector<double>>(natoms, vector<double>(3, 0));
}

void pfHome::pfForce::pfMEAMC::meam_dens_init(Config& cc) {
  getscreen(cc);
  calc_rho1(cc);
}

void pfHome::pfForce::pfMEAMC::getscreen(Config& cc) {
  double sij, rnorm, fc, dfc;
  vector<double> rjk(3, 0);
  double rjk2;
  double xik, xjk, cikj, fcij, sfcij, dfcij, sikj, dfikj;
  double a, Cmax, Cmin, delc, rbound, dCikj;
  double coef1, coef2;
  int elti, eltj, eltk;

  double drinv = 1. / delr_meam;
  for (pfAtom& atm : cc.atoms) {
    elti = atm.tp;

    for (int jj : atm.neighidxHalf) {
      Neigh& ngbj = atm.neighsFull[jj];
      int j = ngbj.aid;
      eltj = cc.atoms[j].tp;

      if (eltj < 0) continue;

      //     First compute screening function itself, sij
      if (ngbj.r > rc_meam) {
        fcij = 0.0;
        dfcij = 0.0;
        sij = 0.0;
      } else {
        rnorm = (rc_meam - ngbj.r) * drinv;
        sij = 1.0;

        //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
        const double rbound = ebound_meam[elti][eltj] * ngbj.r2;
        for (int kk = 0; kk < atm.nneighsFull; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];
          int k = ngbk.aid;
          eltk = cc.atoms[k].tp;
          if (eltk < 0) continue;
          if (k == j) continue;

          for (int it : {0, 1, 2}) rjk[it] = ngbk.dist[it] - ngbj.dist[it];
          rjk2 = square33(rjk);
          if (rjk2 > rbound) continue;
          if (ngbk.r2 > rbound) continue;

          xik = ngbk.r2 / ngbj.r2;
          xjk = rjk2 / ngbj.r2;
          a = 1 - (xik - xjk) * (xik - xjk);
          //     if a < 0, then ellipse equation doesn't describe this case and
          //     atom k can't possibly screen i-j
          if (a <= 0.0) continue;

          cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
          Cmax = Cmax_meam[elti][eltj][eltk];
          Cmin = Cmin_meam[elti][eltj][eltk];
          if (cikj >= Cmax) continue;
          //     note that cikj may be slightly negative (within numerical
          //     tolerance) if atoms are colinear, so don't reject that case
          //     here (other negative cikj cases were handled by the test on "a"
          //     above)
          else if (cikj <= Cmin) {
            sij = 0.0;
            break;
          } else {
            delc = Cmax - Cmin;
            cikj = (cikj - Cmin) / delc;
            sikj = fcut(cikj);
          }
          sij *= sikj;
        }

        fc = dfcut(rnorm, dfc);
        fcij = fc;
        dfcij = dfc * drinv;
      }

      //     Now compute derivatives
      ngbj.dscrfcn = 0.0;
      sfcij = sij * fcij;
      if (iszero(sfcij) || iszero(sfcij - 1.0)) goto LABEL_100;

      rbound = ebound_meam[elti][eltj] * ngbj.r2;

      for (int k = 0; k < atm.nneighsFull; k++) {
        Neigh& ngbk = atm.neighsFull[k];
        eltk = cc.atoms[ngbk.aid].tp;
        if (eltk < 0) continue;
        if (ngbk.aid == ngbj.aid) continue;

        for (int it : {0, 1, 2}) rjk[it] = ngbk.dist[it] - ngbj.dist[it];
        rjk2 = square33(rjk);
        if (rjk2 > rbound) continue;
        if (ngbk.r2 > rbound) continue;

        xik = ngbk.r2 / ngbj.r2;
        xjk = rjk2 / ngbj.r2;
        a = 1 - (xik - xjk) * (xik - xjk);
        //     if a < 0, then ellipse equation doesn't describe this case and
        //     atom k can't possibly screen i-j
        if (a <= 0.0) continue;

        cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
        Cmax = Cmax_meam[elti][eltj][eltk];
        Cmin = Cmin_meam[elti][eltj][eltk];
        if (cikj >= Cmax) {
          continue;
          //     Note that cikj may be slightly negative (within numerical
          //     tolerance) if atoms are colinear, so don't reject that case
          //     here
          //     (other negative cikj cases were handled by the test on "a"
          //     above)
          //     Note that we never have 0<cikj<Cmin here, else sij=0
          //     (rejected above)
        } else {
          delc = Cmax - Cmin;
          cikj = (cikj - Cmin) / delc;
          sikj = dfcut(cikj, dfikj);
          coef1 = dfikj / (delc * sikj);
          dCikj = dCfunc(ngbj.r2, ngbk.r2, rjk2);
          ngbj.dscrfcn = ngbj.dscrfcn + coef1 * dCikj;
        }
      }
      coef1 = sfcij;
      coef2 = sij * dfcij / ngbj.r;
      ngbj.dscrfcn = ngbj.dscrfcn * coef1 - coef2;

    LABEL_100:
      ngbj.scrfcn = sij;
      ngbj.fcpair = fcij;
    }
  }
}

void pfHome::pfForce::pfMEAMC::calc_rho1(Config& cc) {
  int m, n, p, elti, eltj, nv2, nv3;
  double sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;

  for (int ii = 0; ii < cc.natoms; ii++) {
    pfAtom& atm = cc.atoms[ii];
    int i = atm.id;
    elti = atm.tp;

    for (int jj : atm.neighidxHalf) {
      Neigh& ngbj = atm.neighsFull[jj];
      int j = ngbj.aid;

      if (!iszero(ngbj.scrfcn)) {
        sij = ngbj.scrfcn * ngbj.fcpair;

        if (ngbj.r2 < cutforcesq) {
          eltj = cc.atoms[ngbj.aid].tp;

          ai = ngbj.r / re_meam[elti][elti] - 1.0;
          aj = ngbj.r / re_meam[eltj][eltj] - 1.0;
          ro0i = rho0_meam[elti];
          ro0j = rho0_meam[eltj];

          rhoa0j = ro0j * MathSpecial::fm_exp(-beta0_meam[eltj] * aj) * sij;
          rhoa1j = ro0j * MathSpecial::fm_exp(-beta1_meam[eltj] * aj) * sij;
          rhoa2j = ro0j * MathSpecial::fm_exp(-beta2_meam[eltj] * aj) * sij;
          rhoa3j = ro0j * MathSpecial::fm_exp(-beta3_meam[eltj] * aj) * sij;
          rhoa0i = ro0i * MathSpecial::fm_exp(-beta0_meam[elti] * ai) * sij;
          rhoa1i = ro0i * MathSpecial::fm_exp(-beta1_meam[elti] * ai) * sij;
          rhoa2i = ro0i * MathSpecial::fm_exp(-beta2_meam[elti] * ai) * sij;
          rhoa3i = ro0i * MathSpecial::fm_exp(-beta3_meam[elti] * ai) * sij;

          if (ialloy == 1) {
            rhoa1j = rhoa1j * t1_meam[eltj];
            rhoa2j = rhoa2j * t2_meam[eltj];
            rhoa3j = rhoa3j * t3_meam[eltj];
            rhoa1i = rhoa1i * t1_meam[elti];
            rhoa2i = rhoa2i * t2_meam[elti];
            rhoa3i = rhoa3i * t3_meam[elti];
          }
          rho0[i] = rho0[i] + rhoa0j;
          rho0[j] = rho0[j] + rhoa0i;
          // For ialloy = 2, use single-element value (not average)
          if (ialloy != 2) {
            t_ave[i][0] = t_ave[i][0] + t1_meam[eltj] * rhoa0j;
            t_ave[i][1] = t_ave[i][1] + t2_meam[eltj] * rhoa0j;
            t_ave[i][2] = t_ave[i][2] + t3_meam[eltj] * rhoa0j;
            t_ave[j][0] = t_ave[j][0] + t1_meam[elti] * rhoa0i;
            t_ave[j][1] = t_ave[j][1] + t2_meam[elti] * rhoa0i;
            t_ave[j][2] = t_ave[j][2] + t3_meam[elti] * rhoa0i;
          }
          if (ialloy == 1) {
            tsq_ave[i][0] =
                tsq_ave[i][0] + t1_meam[eltj] * t1_meam[eltj] * rhoa0j;
            tsq_ave[i][1] =
                tsq_ave[i][1] + t2_meam[eltj] * t2_meam[eltj] * rhoa0j;
            tsq_ave[i][2] =
                tsq_ave[i][2] + t3_meam[eltj] * t3_meam[eltj] * rhoa0j;
            tsq_ave[j][0] =
                tsq_ave[j][0] + t1_meam[elti] * t1_meam[elti] * rhoa0i;
            tsq_ave[j][1] =
                tsq_ave[j][1] + t2_meam[elti] * t2_meam[elti] * rhoa0i;
            tsq_ave[j][2] =
                tsq_ave[j][2] + t3_meam[elti] * t3_meam[elti] * rhoa0i;
          }
          arho2b[i] = arho2b[i] + rhoa2j;
          arho2b[j] = arho2b[j] + rhoa2i;

          A1j = rhoa1j / ngbj.r;
          A2j = rhoa2j / ngbj.r2;
          A3j = rhoa3j / (ngbj.r2 * ngbj.r);
          A1i = rhoa1i / ngbj.r;
          A2i = rhoa2i / ngbj.r2;
          A3i = rhoa3i / (ngbj.r2 * ngbj.r);
          nv2 = 0;
          nv3 = 0;
          for (m = 0; m < 3; m++) {
            arho1[i][m] = arho1[i][m] + A1j * ngbj.dist[m];
            arho1[j][m] = arho1[j][m] - A1i * ngbj.dist[m];
            arho3b[i][m] = arho3b[i][m] + rhoa3j * ngbj.dist[m] / ngbj.r;
            arho3b[j][m] = arho3b[j][m] - rhoa3i * ngbj.dist[m] / ngbj.r;
            for (n = m; n < 3; n++) {
              arho2[i][nv2] = arho2[i][nv2] + A2j * ngbj.dist[m] * ngbj.dist[n];
              arho2[j][nv2] = arho2[j][nv2] + A2i * ngbj.dist[m] * ngbj.dist[n];
              nv2 = nv2 + 1;
              for (p = n; p < 3; p++) {
                arho3[i][nv3] = arho3[i][nv3] + A3j * ngbj.dist[m] *
                                                    ngbj.dist[n] * ngbj.dist[p];
                arho3[j][nv3] = arho3[j][nv3] - A3i * ngbj.dist[m] *
                                                    ngbj.dist[n] * ngbj.dist[p];
                nv3 = nv3 + 1;
              }
            }
          }
        }
      }
    }
  }
}

void pfHome::pfForce::pfMEAMC::meam_dens_final(Config& cc) {
  int elti, m;
  double rhob, G, dG, Gbar, dGbar, gam, shp[3], Z;
  double B, denom, rho_bkgd;

  //     Complete the calculation of density
  for (int ii = 0; ii < cc.natoms; ii++) {
    pfAtom& atm = cc.atoms[ii];
    int i = atm.id;
    elti = atm.tp;
    if (elti >= 0) {
      rho1[i] = 0.0;
      rho2[i] = -1.0 / 3.0 * arho2b[i] * arho2b[i];
      rho3[i] = 0.0;
      for (m = 0; m < 3; m++) {
        rho1[i] = rho1[i] + arho1[i][m] * arho1[i][m];
        rho3[i] = rho3[i] - 3.0 / 5.0 * arho3b[i][m] * arho3b[i][m];
      }
      for (m = 0; m < 6; m++) {
        rho2[i] = rho2[i] + v2D[m] * arho2[i][m] * arho2[i][m];
      }
      for (m = 0; m < 10; m++) {
        rho3[i] = rho3[i] + v3D[m] * arho3[i][m] * arho3[i][m];
      }

      if (rho0[i] > 0.0) {
        if (ialloy == 1) {
          t_ave[i][0] = t_ave[i][0] / tsq_ave[i][0];
          t_ave[i][1] = t_ave[i][1] / tsq_ave[i][1];
          t_ave[i][2] = t_ave[i][2] / tsq_ave[i][2];
        } else if (ialloy == 2) {
          t_ave[i][0] = t1_meam[elti];
          t_ave[i][1] = t2_meam[elti];
          t_ave[i][2] = t3_meam[elti];
        } else {
          t_ave[i][0] = t_ave[i][0] / rho0[i];
          t_ave[i][1] = t_ave[i][1] / rho0[i];
          t_ave[i][2] = t_ave[i][2] / rho0[i];
        }
      }

      gamma[i] =
          t_ave[i][0] * rho1[i] + t_ave[i][1] * rho2[i] + t_ave[i][2] * rho3[i];

      if (rho0[i] > 0.0) gamma[i] = gamma[i] / (rho0[i] * rho0[i]);

      Z = Z_meam[elti];

      G = G_gam(gamma[i], ibar_meam[elti], errorflag);
      if (errorflag != 0) return;
      get_shpfcn(lattce_meam[elti][elti], shp);
      if (ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
      } else {
        if (mix_ref_t == 1) {
          gam = (t_ave[i][0] * shp[0] + t_ave[i][1] * shp[1] +
                 t_ave[i][2] * shp[2]) /
                (Z * Z);
        } else {
          gam = (t1_meam[elti] * shp[0] + t2_meam[elti] * shp[1] +
                 t3_meam[elti] * shp[2]) /
                (Z * Z);
        }
        Gbar = G_gam(gam, ibar_meam[elti], errorflag);
      }
      rho[i] = rho0[i] * G;

      if (mix_ref_t == 1) {
        if (ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          gam = (t_ave[i][0] * shp[0] + t_ave[i][1] * shp[1] +
                 t_ave[i][2] * shp[2]) /
                (Z * Z);
          Gbar = dG_gam(gam, ibar_meam[elti], dGbar);
        }
        rho_bkgd = rho0_meam[elti] * Z * Gbar;
      } else {
        if (bkgd_dyn == 1) {
          rho_bkgd = rho0_meam[elti] * Z;
        } else {
          rho_bkgd = rho_ref_meam[elti];
        }
      }
      rhob = rho[i] / rho_bkgd;
      denom = 1.0 / rho_bkgd;

      G = dG_gam(gamma[i], ibar_meam[elti], dG);

      dgamma1[i] = (G - 2 * dG * gamma[i]) * denom;

      if (!iszero(rho0[i])) {
        dgamma2[i] = (dG / rho0[i]) * denom;
      } else {
        dgamma2[i] = 0.0;
      }

      //     dgamma3 is nonzero only if we are using the "mixed" rule for
      //     computing t in the reference system (which is not correct, but
      //     included for backward compatibility
      if (mix_ref_t == 1) {
        dgamma3[i] = rho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
      } else {
        dgamma3[i] = 0.0;
      }

      B = A_meam[elti] * Ec_meam[elti][elti];

      if (!iszero(rhob)) {
        if (emb_lin_neg == 1 && rhob <= 0) {
          frhop[i] = -B;
          cc.fitengy -= B * rhob;
          atm.eng -= B * rhob;
        } else {
          frhop[i] = B * (log(rhob) + 1.0);
          cc.fitengy += B * rhob * log(rhob);
          atm.eng += B * rhob * log(rhob);
        }

      } else {
        if (emb_lin_neg == 1) {
          frhop[i] = -B;
        } else {
          frhop[i] = B;
        }
      }
    }
  }
}
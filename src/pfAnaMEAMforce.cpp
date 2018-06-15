/*
 * @Author: chaomy
 * @Date:   2018-02-02 15:13:57
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 00:17:15
 */

#include "pfMEAMC.h"

void pfHome::pfForce::pfMEAMC::meam_force(Config& cc) {
  int kk, m, n, p, q;
  int nv2, nv3, elti, eltj, eltk, ind;
  double delij[3], rij2, rij, rij3;
  double v[6], fi[3], fj[3];
  double third, sixth;
  double pp, dUdrij, dUdsij, dUdrijm[3], force, forcem;
  double r, recip, phi, phip;
  double sij;
  double a1, a1i, a1j, a2, a2i, a2j;
  double a3i, a3j;
  double shpi[3], shpj[3];
  double ai, aj, ro0i, ro0j, invrei, invrej;
  double rhoa0j, drhoa0j, rhoa0i, drhoa0i;
  double rhoa1j, drhoa1j, rhoa1i, drhoa1i;
  double rhoa2j, drhoa2j, rhoa2i, drhoa2i;
  double a3, a3a, rhoa3j, drhoa3j, rhoa3i, drhoa3i;
  double drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  double drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  double drho1drm1[3], drho1drm2[3];
  double drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  double drho2drm1[3], drho2drm2[3];
  double drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  double drho3drm1[3], drho3drm2[3];
  double dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  double dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  double dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  double drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3], drhodrm2[3];
  double arg;
  double arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  double dsij1, dsij2, force1, force2;
  double t1i, t2i, t3i, t1j, t2j, t3j;

  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  //     Compute forces atom i
  for (pfAtom& atm : cc.atoms) {
    int i = atm.id;
    elti = atm.tp;
    if (elti < 0) return;

    for (int jn : atm.neighidxHalf) {
      Neigh& ngbj = atm.neighsFull[jn];
      int j = ngbj.aid;
      eltj = cc.atoms[j].tp;

      if (!iszero(ngbj.scrfcn) && (eltj >= 0)) {
        sij = ngbj.scrfcn * ngbj.fcpair;
        delij[0] = ngbj.dist[0];
        delij[1] = ngbj.dist[1];
        delij[2] = ngbj.dist[2];
        rij2 = ngbj.r2;
        if (rij2 < cutforcesq) {
          r = rij = ngbj.r;

          //     Compute phi and phip
          ind = eltind[elti][eltj];
          pp = rij * rdrar;
          kk = (int)pp;
          kk = std::min(kk, nrar - 2);
          pp = pp - kk;
          pp = std::min(pp, 1.0);
          phi = ((phirar3[ind][kk] * pp + phirar2[ind][kk]) * pp +
                 phirar1[ind][kk]) *
                    pp +
                phirar[ind][kk];
          phip = (phirar6[ind][kk] * pp + phirar5[ind][kk]) * pp +
                 phirar4[ind][kk];
          recip = 1.0 / r;

          cc.fitengy += phi * sij;
          atm.eng += 0.5 * phi * sij;
          cc.atoms[j].eng += 0.5 * phi * sij;

          //     Compute pair densities and derivatives
          invrei = 1.0 / re_meam[elti][elti];
          ai = rij * invrei - 1.0;
          ro0i = rho0_meam[elti];
          rhoa0i = ro0i * MathSpecial::fm_exp(-beta0_meam[elti] * ai);
          drhoa0i = -beta0_meam[elti] * invrei * rhoa0i;
          rhoa1i = ro0i * MathSpecial::fm_exp(-beta1_meam[elti] * ai);
          drhoa1i = -beta1_meam[elti] * invrei * rhoa1i;
          rhoa2i = ro0i * MathSpecial::fm_exp(-beta2_meam[elti] * ai);
          drhoa2i = -beta2_meam[elti] * invrei * rhoa2i;
          rhoa3i = ro0i * MathSpecial::fm_exp(-beta3_meam[elti] * ai);
          drhoa3i = -beta3_meam[elti] * invrei * rhoa3i;

          if (elti != eltj) {
            invrej = 1.0 / re_meam[eltj][eltj];
            aj = rij * invrej - 1.0;
            ro0j = rho0_meam[eltj];
            rhoa0j = ro0j * MathSpecial::fm_exp(-beta0_meam[eltj] * aj);
            drhoa0j = -beta0_meam[eltj] * invrej * rhoa0j;
            rhoa1j = ro0j * MathSpecial::fm_exp(-beta1_meam[eltj] * aj);
            drhoa1j = -beta1_meam[eltj] * invrej * rhoa1j;
            rhoa2j = ro0j * MathSpecial::fm_exp(-beta2_meam[eltj] * aj);
            drhoa2j = -beta2_meam[eltj] * invrej * rhoa2j;
            rhoa3j = ro0j * MathSpecial::fm_exp(-beta3_meam[eltj] * aj);
            drhoa3j = -beta3_meam[eltj] * invrej * rhoa3j;
          } else {
            rhoa0j = rhoa0i;
            drhoa0j = drhoa0i;
            rhoa1j = rhoa1i;
            drhoa1j = drhoa1i;
            rhoa2j = rhoa2i;
            drhoa2j = drhoa2i;
            rhoa3j = rhoa3i;
            drhoa3j = drhoa3i;
          }

          const double t1mi = t1_meam[elti];
          const double t2mi = t2_meam[elti];
          const double t3mi = t3_meam[elti];
          const double t1mj = t1_meam[eltj];
          const double t2mj = t2_meam[eltj];
          const double t3mj = t3_meam[eltj];

          if (ialloy == 1) {
            rhoa1j *= t1mj;
            rhoa2j *= t2mj;
            rhoa3j *= t3mj;
            rhoa1i *= t1mi;
            rhoa2i *= t2mi;
            rhoa3i *= t3mi;
            drhoa1j *= t1mj;
            drhoa2j *= t2mj;
            drhoa3j *= t3mj;
            drhoa1i *= t1mi;
            drhoa2i *= t2mi;
            drhoa3i *= t3mi;
          }

          nv2 = 0;
          nv3 = 0;
          arg1i1 = 0.0;
          arg1j1 = 0.0;
          arg1i2 = 0.0;
          arg1j2 = 0.0;
          arg1i3 = 0.0;
          arg1j3 = 0.0;
          arg3i3 = 0.0;
          arg3j3 = 0.0;
          for (n = 0; n < 3; n++) {
            for (p = n; p < 3; p++) {
              for (q = p; q < 3; q++) {
                arg = delij[n] * delij[p] * delij[q] * v3D[nv3];
                arg1i3 = arg1i3 + arho3[i][nv3] * arg;
                arg1j3 = arg1j3 - arho3[j][nv3] * arg;
                nv3 = nv3 + 1;
              }
              arg = delij[n] * delij[p] * v2D[nv2];
              arg1i2 = arg1i2 + arho2[i][nv2] * arg;
              arg1j2 = arg1j2 + arho2[j][nv2] * arg;
              nv2 = nv2 + 1;
            }
            arg1i1 = arg1i1 + arho1[i][n] * delij[n];
            arg1j1 = arg1j1 - arho1[j][n] * delij[n];
            arg3i3 = arg3i3 + arho3b[i][n] * delij[n];
            arg3j3 = arg3j3 - arho3b[j][n] * delij[n];
          }

          //     rho0 terms
          drho0dr1 = drhoa0j * sij;
          drho0dr2 = drhoa0i * sij;

          //     rho1 terms
          a1 = 2 * sij / rij;
          drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1;
          drho1dr2 = a1 * (drhoa1i - rhoa1i / rij) * arg1j1;
          a1 = 2.0 * sij / rij;
          for (m = 0; m < 3; m++) {
            drho1drm1[m] = a1 * rhoa1j * arho1[i][m];
            drho1drm2[m] = -a1 * rhoa1i * arho1[j][m];
          }

          //     rho2 terms
          a2 = 2 * sij / rij2;
          drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 -
                     2.0 / 3.0 * arho2b[i] * drhoa2j * sij;
          drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 -
                     2.0 / 3.0 * arho2b[j] * drhoa2i * sij;
          a2 = 4 * sij / rij2;
          for (m = 0; m < 3; m++) {
            drho2drm1[m] = 0.0;
            drho2drm2[m] = 0.0;
            for (n = 0; n < 3; n++) {
              drho2drm1[m] = drho2drm1[m] + arho2[i][vind2D[m][n]] * delij[n];
              drho2drm2[m] = drho2drm2[m] - arho2[j][vind2D[m][n]] * delij[n];
            }
            drho2drm1[m] = a2 * rhoa2j * drho2drm1[m];
            drho2drm2[m] = -a2 * rhoa2i * drho2drm2[m];
          }

          //     rho3 terms
          rij3 = rij * rij2;
          a3 = 2 * sij / rij3;
          a3a = 6.0 / 5.0 * sij / rij;
          drho3dr1 = a3 * (drhoa3j - 3 * rhoa3j / rij) * arg1i3 -
                     a3a * (drhoa3j - rhoa3j / rij) * arg3i3;
          drho3dr2 = a3 * (drhoa3i - 3 * rhoa3i / rij) * arg1j3 -
                     a3a * (drhoa3i - rhoa3i / rij) * arg3j3;
          a3 = 6 * sij / rij3;
          a3a = 6 * sij / (5 * rij);
          for (m = 0; m < 3; m++) {
            drho3drm1[m] = 0.0;
            drho3drm2[m] = 0.0;
            nv2 = 0;
            for (n = 0; n < 3; n++) {
              for (p = n; p < 3; p++) {
                arg = delij[n] * delij[p] * v2D[nv2];
                drho3drm1[m] = drho3drm1[m] + arho3[i][vind3D[m][n][p]] * arg;
                drho3drm2[m] = drho3drm2[m] + arho3[j][vind3D[m][n][p]] * arg;
                nv2 = nv2 + 1;
              }
            }
            drho3drm1[m] = (a3 * drho3drm1[m] - a3a * arho3b[i][m]) * rhoa3j;
            drho3drm2[m] = (-a3 * drho3drm2[m] + a3a * arho3b[j][m]) * rhoa3i;
          }

          //     Compute derivatives of weighting functions t wrt rij
          t1i = t_ave[i][0];
          t2i = t_ave[i][1];
          t3i = t_ave[i][2];
          t1j = t_ave[j][0];
          t2j = t_ave[j][1];
          t3j = t_ave[j][2];

          if (ialloy == 1) {
            a1i = 0.0;
            a1j = 0.0;
            a2i = 0.0;
            a2j = 0.0;
            a3i = 0.0;
            a3j = 0.0;
            if (!iszero(tsq_ave[i][0])) a1i = drhoa0j * sij / tsq_ave[i][0];
            if (!iszero(tsq_ave[j][0])) a1j = drhoa0i * sij / tsq_ave[j][0];
            if (!iszero(tsq_ave[i][1])) a2i = drhoa0j * sij / tsq_ave[i][1];
            if (!iszero(tsq_ave[j][1])) a2j = drhoa0i * sij / tsq_ave[j][1];
            if (!iszero(tsq_ave[i][2])) a3i = drhoa0j * sij / tsq_ave[i][2];
            if (!iszero(tsq_ave[j][2])) a3j = drhoa0i * sij / tsq_ave[j][2];

            dt1dr1 = a1i * (t1mj - t1i * t1mj * t1mj);
            dt1dr2 = a1j * (t1mi - t1j * t1mi * t1mi);
            dt2dr1 = a2i * (t2mj - t2i * t2mj * t2mj);
            dt2dr2 = a2j * (t2mi - t2j * t2mi * t2mi);
            dt3dr1 = a3i * (t3mj - t3i * t3mj * t3mj);
            dt3dr2 = a3j * (t3mi - t3j * t3mi * t3mi);

          } else if (ialloy == 2) {
            dt1dr1 = 0.0;
            dt1dr2 = 0.0;
            dt2dr1 = 0.0;
            dt2dr2 = 0.0;
            dt3dr1 = 0.0;
            dt3dr2 = 0.0;

          } else {
            ai = 0.0;
            if (!iszero(rho0[i])) ai = drhoa0j * sij / rho0[i];
            aj = 0.0;
            if (!iszero(rho0[j])) aj = drhoa0i * sij / rho0[j];

            dt1dr1 = ai * (t1mj - t1i);
            dt1dr2 = aj * (t1mi - t1j);
            dt2dr1 = ai * (t2mj - t2i);
            dt2dr2 = aj * (t2mi - t2j);
            dt3dr1 = ai * (t3mj - t3i);
            dt3dr2 = aj * (t3mi - t3j);
          }

          //     Compute derivatives of total density wrt rij, sij and rij(3)
          get_shpfcn(lattce_meam[elti][elti], shpi);
          get_shpfcn(lattce_meam[eltj][eltj], shpj);
          drhodr1 = dgamma1[i] * drho0dr1 +
                    dgamma2[i] *
                        (dt1dr1 * rho1[i] + t1i * drho1dr1 + dt2dr1 * rho2[i] +
                         t2i * drho2dr1 + dt3dr1 * rho3[i] + t3i * drho3dr1) -
                    dgamma3[i] * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 +
                                  shpi[2] * dt3dr1);
          drhodr2 = dgamma1[j] * drho0dr2 +
                    dgamma2[j] *
                        (dt1dr2 * rho1[j] + t1j * drho1dr2 + dt2dr2 * rho2[j] +
                         t2j * drho2dr2 + dt3dr2 * rho3[j] + t3j * drho3dr2) -
                    dgamma3[j] * (shpj[0] * dt1dr2 + shpj[1] * dt2dr2 +
                                  shpj[2] * dt3dr2);
          for (m = 0; m < 3; m++) {
            drhodrm1[m] = 0.0;
            drhodrm2[m] = 0.0;
            drhodrm1[m] =
                dgamma2[i] *
                (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]);
            drhodrm2[m] =
                dgamma2[j] *
                (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
          }

          //     Compute derivatives wrt sij, but only if necessary
          if (!iszero(ngbj.dscrfcn)) {
            drho0ds1 = rhoa0j;
            drho0ds2 = rhoa0i;
            a1 = 2.0 / rij;
            drho1ds1 = a1 * rhoa1j * arg1i1;
            drho1ds2 = a1 * rhoa1i * arg1j1;
            a2 = 2.0 / rij2;
            drho2ds1 = a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * arho2b[i] * rhoa2j;
            drho2ds2 = a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * arho2b[j] * rhoa2i;
            a3 = 2.0 / rij3;
            a3a = 6.0 / (5.0 * rij);
            drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
            drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

            if (ialloy == 1) {
              a1i = 0.0;
              a1j = 0.0;
              a2i = 0.0;
              a2j = 0.0;
              a3i = 0.0;
              a3j = 0.0;
              if (!iszero(tsq_ave[i][0])) a1i = rhoa0j / tsq_ave[i][0];
              if (!iszero(tsq_ave[j][0])) a1j = rhoa0i / tsq_ave[j][0];
              if (!iszero(tsq_ave[i][1])) a2i = rhoa0j / tsq_ave[i][1];
              if (!iszero(tsq_ave[j][1])) a2j = rhoa0i / tsq_ave[j][1];
              if (!iszero(tsq_ave[i][2])) a3i = rhoa0j / tsq_ave[i][2];
              if (!iszero(tsq_ave[j][2])) a3j = rhoa0i / tsq_ave[j][2];

              dt1ds1 = a1i * (t1mj - t1i * pow(t1mj, 2));
              dt1ds2 = a1j * (t1mi - t1j * pow(t1mi, 2));
              dt2ds1 = a2i * (t2mj - t2i * pow(t2mj, 2));
              dt2ds2 = a2j * (t2mi - t2j * pow(t2mi, 2));
              dt3ds1 = a3i * (t3mj - t3i * pow(t3mj, 2));
              dt3ds2 = a3j * (t3mi - t3j * pow(t3mi, 2));

            } else if (ialloy == 2) {
              dt1ds1 = 0.0;
              dt1ds2 = 0.0;
              dt2ds1 = 0.0;
              dt2ds2 = 0.0;
              dt3ds1 = 0.0;
              dt3ds2 = 0.0;

            } else {
              ai = 0.0;
              if (!iszero(rho0[i])) ai = rhoa0j / rho0[i];
              aj = 0.0;
              if (!iszero(rho0[j])) aj = rhoa0i / rho0[j];

              dt1ds1 = ai * (t1mj - t1i);
              dt1ds2 = aj * (t1mi - t1j);
              dt2ds1 = ai * (t2mj - t2i);
              dt2ds2 = aj * (t2mi - t2j);
              dt3ds1 = ai * (t3mj - t3i);
              dt3ds2 = aj * (t3mi - t3j);
            }

            drhods1 = dgamma1[i] * drho0ds1 +
                      dgamma2[i] * (dt1ds1 * rho1[i] + t1i * drho1ds1 +
                                    dt2ds1 * rho2[i] + t2i * drho2ds1 +
                                    dt3ds1 * rho3[i] + t3i * drho3ds1) -
                      dgamma3[i] * (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 +
                                    shpi[2] * dt3ds1);
            drhods2 = dgamma1[j] * drho0ds2 +
                      dgamma2[j] * (dt1ds2 * rho1[j] + t1j * drho1ds2 +
                                    dt2ds2 * rho2[j] + t2j * drho2ds2 +
                                    dt3ds2 * rho3[j] + t3j * drho3ds2) -
                      dgamma3[j] * (shpj[0] * dt1ds2 + shpj[1] * dt2ds2 +
                                    shpj[2] * dt3ds2);
          }

          //     Compute derivatives of energy wrt rij, sij and rij[3]
          dUdrij = phip * sij + frhop[i] * drhodr1 + frhop[j] * drhodr2;
          dUdsij = 0.0;
          if (!iszero(ngbj.dscrfcn)) {
            dUdsij = phi + frhop[i] * drhods1 + frhop[j] * drhods2;
          }
          for (m = 0; m < 3; m++) {
            dUdrijm[m] = frhop[i] * drhodrm1[m] + frhop[j] * drhodrm2[m];
          }

          //     Add the part of the force due to dUdrij and dUdsij

          force = dUdrij * recip + dUdsij * ngbj.dscrfcn;
          for (m = 0; m < 3; m++) {
            forcem = delij[m] * force + dUdrijm[m];
            atm.fitfrc[m] += forcem;
            cc.atoms[j].fitfrc[m] -= forcem;
          }

          //     Tabulate per-atom virial as symmetrized stress tensor

          fi[0] = delij[0] * force + dUdrijm[0];
          fi[1] = delij[1] * force + dUdrijm[1];
          fi[2] = delij[2] * force + dUdrijm[2];
          v[0] = -0.5 * (delij[0] * fi[0]);
          v[1] = -0.5 * (delij[1] * fi[1]);
          v[2] = -0.5 * (delij[2] * fi[2]);
          v[3] = -0.25 * (delij[0] * fi[1] + delij[1] * fi[0]);
          v[4] = -0.25 * (delij[0] * fi[2] + delij[2] * fi[0]);
          v[5] = -0.25 * (delij[1] * fi[2] + delij[2] * fi[1]);

          for (m = 0; m < 6; m++) {
            atm.sts[m] += v[m];
            cc.atoms[j].sts[m] += v[m];
          }

          //     Now compute forces on other atoms k due to change in sij

          if (iszero(sij) || iszero(sij - 1.0)) continue;  //: cont jn loop

          double dxik(0), dyik(0), dzik(0);
          double dxjk(0), dyjk(0), dzjk(0);

          for (int kn = 0; kn < atm.nneighsFull; kn++) {
            Neigh& ngbk = atm.neighsFull[kn];
            int k = ngbk.aid;
            eltk = cc.atoms[k].tp;

            if (k != j && eltk >= 0) {
              double xik, xjk, cikj, sikj, dfc, a;
              double dCikj1, dCikj2;
              double delc, rik2, rjk2;

              sij = ngbj.scrfcn * ngbj.fcpair;
              const double Cmax = Cmax_meam[elti][eltj][eltk];
              const double Cmin = Cmin_meam[elti][eltj][eltk];

              dsij1 = 0.0;
              dsij2 = 0.0;
              if (!iszero(sij) && !iszero(sij - 1.0)) {
                const double rbound = rij2 * ebound_meam[elti][eltj];
                delc = Cmax - Cmin;
                dxjk = ngbk.dist[0] - ngbj.dist[0];
                dyjk = ngbk.dist[1] - ngbj.dist[1];
                dzjk = ngbk.dist[2] - ngbj.dist[2];
                rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
                if (rjk2 <= rbound) {
                  rik2 = ngbk.r2;
                  if (rik2 <= rbound) {
                    dxik = ngbk.dist[0];
                    dyik = ngbk.dist[1];
                    dzik = ngbk.dist[2];
                    xik = rik2 / rij2;
                    xjk = rjk2 / rij2;
                    a = 1 - (xik - xjk) * (xik - xjk);
                    if (!iszero(a)) {
                      cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
                      if (cikj >= Cmin && cikj <= Cmax) {
                        cikj = (cikj - Cmin) / delc;
                        sikj = dfcut(cikj, dfc);
                        dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
                        a = sij / delc * dfc / sikj;
                        dsij1 = a * dCikj1;
                        dsij2 = a * dCikj2;
                      }
                    }
                  }
                }
              }

              if (!iszero(dsij1) || !iszero(dsij2)) {
                force1 = dUdsij * dsij1;
                force2 = dUdsij * dsij2;

                atm.fitfrc[0] += force1 * dxik;
                atm.fitfrc[1] += force1 * dyik;
                atm.fitfrc[2] += force1 * dzik;
                cc.atoms[j].fitfrc[0] += force2 * dxjk;
                cc.atoms[j].fitfrc[1] += force2 * dyjk;
                cc.atoms[j].fitfrc[2] += force2 * dzjk;
                cc.atoms[k].fitfrc[0] -= force1 * dxik + force2 * dxjk;
                cc.atoms[k].fitfrc[1] -= force1 * dyik + force2 * dyjk;
                cc.atoms[k].fitfrc[2] -= force1 * dzik + force2 * dzjk;

                //     Tabulate per-atom virial as symmetrized stress tensor

                fi[0] = force1 * dxik;
                fi[1] = force1 * dyik;
                fi[2] = force1 * dzik;
                fj[0] = force2 * dxjk;
                fj[1] = force2 * dyjk;
                fj[2] = force2 * dzjk;
                v[0] = -third * (dxik * fi[0] + dxjk * fj[0]);
                v[1] = -third * (dyik * fi[1] + dyjk * fj[1]);
                v[2] = -third * (dzik * fi[2] + dzjk * fj[2]);
                v[3] = -sixth * (dxik * fi[1] + dxjk * fj[1] + dyik * fi[0] +
                                 dyjk * fj[0]);
                v[4] = -sixth * (dxik * fi[2] + dxjk * fj[2] + dzik * fi[0] +
                                 dzjk * fj[0]);
                v[5] = -sixth * (dyik * fi[2] + dyjk * fj[2] + dzik * fi[1] +
                                 dzjk * fj[1]);

                for (m = 0; m < 6; m++) {
                  atm.sts[m] += v[m];
                  cc.atoms[j].sts[m] += v[m];
                  cc.atoms[k].sts[m] += v[m];
                }
              }
            }  //     end of k loop
          }
        }
      }  //     end of j loop
    }
  }
}

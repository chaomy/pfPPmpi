/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-26 16:03:18
 */

#include "pfForce.h"
#include "pfLmpDrv.h"

double pfHome::pfForce::forceMEAMSStressPunish(const arma::mat &vv, int tg) {
  while (true) {
    broadcast(cmm, tg, PFROOT);
    if (tg == EXT) break;

    int cnt = 0;
    for (int i : {0, 1, 2, 3, 4}) {
      if (optidx[i] == 0) continue;
      Func &ff = funcs[i];
      for (int j : ff.rlxid) ff.yy[j] = vv[cnt++];
    }

    for (int i = 0; i < funcs.size(); i++) {  // broadcast functions
      broadcast(cmm, funcs[i].xx, PFROOT);
      broadcast(cmm, funcs[i].yy, PFROOT);
    }

    for (Func &ff : funcs) ff.s.set_points(ff.xx, ff.yy);

    double efrc = 0.0, eengy = 0.0, estrs = 0.0;
    error.frc = 0.0, error.pnsh = 0.0;
    omaxrho = -1e10, ominrho = 1e10;

    // int ww = 1;
    // for (int it : smthidx) {  // covarance of second derivative
    //   vector<double> &vv = funcs[it].s.m_b;
    //   double mn = 0.0, cov = 0.0;
    //   for (int i = ww + 1; i < vv.size() - ww; i++) {
    //     for (int it = -ww; it <= ww; it++) mn += vv[i + it];
    //     mn /= (2 * ww + 1);
    //     for (int it = -ww; it <= ww; it++) cov += square11(vv[i + it] - mn);
    //   }
    //   error.pnsh += cov / (2 * ww + 1);
    //   // error.pnsh += square11(mn);
    // }

    for (int i : locls) {
      Config &cnf = configs[i];
      forceMEAMS(cnf);
      for (pfAtom &atm : cnf.atoms) {
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it];
          efrc += cnf.weigh *
                  square11((atm.fitfrc[it] - atm.frc[it]) * atm.fweigh[it]);
        }
      }
      omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
      ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
      eengy += cnf.weigh * square11(cnf.fitengy - cnf.engy);

      for (int it : {0, 1, 2, 3, 4, 5})
        estrs += square11(cnf.fitstrs[it] - cnf.strs[it]);
      estrs *= cnf.weigh;
    }

    error.pnsh *= weigh.pnsh;
    estrs *= weigh.stss;
    eengy *= weigh.engy;

    reduce(cmm, eengy, error.engy, std::plus<double>(), PFROOT);
    reduce(cmm, estrs, error.stss, std::plus<double>(), PFROOT);
    reduce(cmm, efrc, error.frc, std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  return error.frc + error.engy + error.stss + error.pnsh;
}

double pfHome::pfForce::forceMEAMSStress(const arma::mat &vv, int tg) {
  while (true) {
    broadcast(cmm, tg, PFROOT);
    if (tg == EXT) break;

    int cnt = 0;
    for (int i : {0, 1, 2, 3, 4}) {
      if (optidx[i] == 0) continue;
      Func &ff = funcs[i];
      for (int j : ff.rlxid) ff.yy[j] = vv[cnt++];
    }

    for (int i = 0; i < funcs.size(); i++) {  // broadcast functions
      broadcast(cmm, funcs[i].xx, PFROOT);
      broadcast(cmm, funcs[i].yy, PFROOT);
    }

    for (Func &ff : funcs) ff.s.set_points(ff.xx, ff.yy);

    double efrc = 0.0, eengy = 0.0, estrs = 0.0;
    error.frc = 0.0;
    omaxrho = -1e10, ominrho = 1e10;

    for (int i : locls) {
      Config &cnf = configs[i];
      forceMEAMS(cnf);
      for (pfAtom &atm : cnf.atoms) {
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it];
          efrc += cnf.weigh *
                  square11((atm.fitfrc[it] - atm.frc[it]) * atm.fweigh[it]);
        }
      }
      eengy += cnf.weigh * square11(cnf.fitengy - cnf.engy);
      omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
      ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;

      for (int it : {0, 1, 2, 3, 4, 5})
        estrs += square11(cnf.fitstrs[it] - cnf.strs[it]);
      estrs *= cnf.weigh;
    }
    reduce(cmm, weigh.engy * eengy, error.engy, std::plus<double>(), PFROOT);
    reduce(cmm, weigh.stss * estrs, error.stss, std::plus<double>(), PFROOT);
    reduce(cmm, efrc, error.frc, std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  return error.frc + error.engy + error.stss;
}

double pfHome::pfForce::forceMEAMS(const arma::mat &vv) {
  error.frc = 0.0, error.pnsh = 0.0;
  omaxrho = -1e10, ominrho = 1e10;

  int cnt = 0;
  for (int i = 0; i < funcs.size(); i++) {
    Func &ff = funcs[i];
    if (i == PHI || i == RHO || i == MEAMF) {
      for (int j = 0; j < ff.npts - 1; j++) ff.yy[j] = vv[cnt++];
    } else {
      for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];
    }
    ff.s.set_points(ff.xx, ff.yy);
  }

  int ls[] = {PHI, RHO, MEAMF};
  for (int it : ls) {
    for (double ee : funcs[it].s.m_b) error.pnsh += square11(ee);
    for (int i = 0; i < funcs[it].s.m_b.size() - 1; i++)
      error.pnsh += 0.5 * funcs[it].s.m_b[i] * funcs[it].s.m_b[i + 1];
  }

  error.frc = 0.0;
  for (Config &cnf : configs) {
    forceMEAMS(cnf);
    for (pfAtom &atm : cnf.atoms)
      for (int it : {X, Y, Z}) {
        atm.fitfrc[it] =
            atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it] - atm.frc[it];
        error.frc += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }

    ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
    omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
    error.frc += square11(cnf.fitengy - cnf.engy);
  }  // cc
  error.pnsh *= weigh.pnsh;
  error.frc *= 1e2;
  return error.frc + error.pnsh;
}

void pfHome::pfForce::forceMEAMSStress(Config &cnf) {
  cnf.phiengy = cnf.emfengy = 0.0;
  cnf.rhomi = 1e10, cnf.rhomx = -1e10;
  for (auto &ee : cnf.fitstrs) ee = 0.0;
  for (pfAtom &atm : cnf.atoms) { /* loop over atoms to reset values */
    atm.crho = 0.0;
    for (int it : {X, Y, Z})
      atm.phifrc[it] = atm.rhofrc[it] = atm.trifrc[it] = 0.0;
  }  // ii
  double e0 = funcs[EMF].s(0.0);
  for (pfAtom &atm : cnf.atoms) { /* loop over atoms pairs, densities */
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];

      // pair (phi)
      funcs[PHI].s.deriv(ngbj.slots[PHI], ngbj.shifts[PHI], ngbj.phi,
                         ngbj.phig);
      cnf.phiengy += ngbj.phi;

      double tmp[3];
      for (int it : {X, Y, Z})
        atm.phifrc[it] += (tmp[it] = ngbj.dist2r[it] * ngbj.phig);

      cnf.fitstrs[XX] -= 0.5 * ngbj.dist[X] * tmp[X];
      cnf.fitstrs[YY] -= 0.5 * ngbj.dist[Y] * tmp[Y];
      cnf.fitstrs[ZZ] -= 0.5 * ngbj.dist[Z] * tmp[Z];
      cnf.fitstrs[XY] -= 0.5 * ngbj.dist[X] * tmp[Y];
      cnf.fitstrs[YZ] -= 0.5 * ngbj.dist[Y] * tmp[Z];
      cnf.fitstrs[ZX] -= 0.5 * ngbj.dist[Z] * tmp[X];

      // rho
      funcs[RHO].s.deriv(ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.rho,
                         ngbj.rhog);
      atm.crho += ngbj.rho;
      // partial sum
      double psum = 0.0;
      funcs[MEAMF].s.deriv(ngbj.slots[2], ngbj.shifts[2], ngbj.fval,
                           ngbj.fgrad);
      for (int kk = 0; kk < jj; kk++) {
        Neigh &ngbk = atm.neighsFull[kk];

        Angle &agl = atm.angMat[jj][kk];
        funcs[MEAMG].s.deriv(agl.slot, agl.shift, agl.gval, agl.ggrad);
        psum += ngbk.fval * agl.gval;
      }  // kk
      atm.crho += psum * ngbj.fval;
    }  // jj

    double embE;
    funcs[EMF].s.deriv(atm.crho, embE, atm.gradF);

    cnf.emfengy += (embE - e0);

    cnf.rhomx = atm.crho > cnf.rhomx ? atm.crho : cnf.rhomx;
    cnf.rhomi = atm.crho < cnf.rhomi ? atm.crho : cnf.rhomi;

    double forces_i[3] = {0.0, 0.0, 0.0};
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];

      double f_rij = ngbj.fval;
      double f_rij_prime = ngbj.fgrad;

      double forces_j[3] = {0.0, 0.0, 0.0};
      for (int kk = 0; kk < jj; kk++) {
        Neigh &ngbk = atm.neighsFull[kk];

        double gcos = atm.angMat[jj][kk].gcos;
        double gval = atm.angMat[jj][kk].gval;
        double gprime = atm.angMat[jj][kk].ggrad;

        double f_rik = ngbk.fval;
        double f_rik_prime = ngbk.fgrad;

        double fij = -atm.gradF * gval * f_rik * f_rij_prime;
        double fik = -atm.gradF * gval * f_rij * f_rik_prime;

        double prefactor = atm.gradF * f_rij * f_rik * gprime;

        double prefactor_ij = prefactor / ngbj.r;
        double prefactor_ik = prefactor / ngbk.r;

        fij += prefactor_ij * gcos;
        fik += prefactor_ik * gcos;

        double fj[3], fk[3];

        for (int it : {X, Y, Z}) {
          fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
          forces_j[it] += fj[it];

          fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
          forces_i[it] -= fk[it];
          cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
        }

        // stresses
        cnf.fitstrs[XX] += ngbj.dist[X] * fj[X] + ngbk.dist[X] * fk[X];
        cnf.fitstrs[YY] += ngbj.dist[Y] * fj[Y] + ngbk.dist[Y] * fk[Y];
        cnf.fitstrs[ZZ] += ngbj.dist[Z] * fj[Z] + ngbk.dist[Z] * fk[Z];
        cnf.fitstrs[XY] += ngbj.dist[X] * fj[Y] + ngbk.dist[X] * fk[Y];
        cnf.fitstrs[YZ] += ngbj.dist[Y] * fj[Z] + ngbk.dist[Y] * fk[Z];
        cnf.fitstrs[ZX] += ngbj.dist[Z] * fj[X] + ngbk.dist[Z] * fk[X];
      }  // loop over kk
      for (int it : {X, Y, Z}) {
        atm.trifrc[it] -= forces_j[it];
        cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
      }
    }  // loop over jj
    for (int it : {X, Y, Z}) atm.trifrc[it] += forces_i[it];
  }                             // atm
  for (pfAtom &atm : cnf.atoms) /* eambedding forces */
    for (Neigh &ngb : atm.neighsFull) {
      // if (ngb.r < funcs[RHO].xx.back()) {
      double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
      double tmp[3];
      for (int it : {X, Y, Z})
        atm.rhofrc[it] += (tmp[it] = ngb.dist2r[it] * emf);
      cnf.fitstrs[XX] -= 0.5 * ngb.dist[X] * tmp[X];
      cnf.fitstrs[YY] -= 0.5 * ngb.dist[Y] * tmp[Y];
      cnf.fitstrs[ZZ] -= 0.5 * ngb.dist[Z] * tmp[Z];
      cnf.fitstrs[XY] -= 0.5 * ngb.dist[X] * tmp[Y];
      cnf.fitstrs[YZ] -= 0.5 * ngb.dist[Y] * tmp[Z];
      cnf.fitstrs[ZX] -= 0.5 * ngb.dist[Z] * tmp[X];
      // }  // cutoff(rho)
    }  // nn
  cnf.fitengy = (cnf.phiengy / 2. + cnf.emfengy) / cnf.natoms;
  for (auto &ee : cnf.fitstrs) ee /= cnf.vol;
}

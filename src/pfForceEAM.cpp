/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-20 14:35:57
 */

#include "pfHome.h"

double pfHome::forceEAM(const arma::mat& vv) {
  error["frc"] = 0.0, error["punish"] = 0.0, error["shift"] = 0.0;
  omaxrho = -1e10, ominrho = 1e10;
  int cnt = 0;

  for (int i = 0; i < nfuncs; i++) { /* interpolates */
    Func& ff = funcs[i];
    double mxf = -1e10, mif = 1e10;
    int nt = (i == PHI || i == RHO) ? ff.npts - 1 : ff.npts;
    for (int j = 0; j < nt; j++) {
      ff.yy[j] = vv[cnt++];
      mxf = ff.yy[j] > mxf ? ff.yy[j] : mxf;
      mif = ff.yy[j] < mif ? ff.yy[j] : mif;
    }
    ff.s.set_points(ff.xx, ff.yy);
    ff.rng = mxf - mif;
  }

  /* reverse */
  // for (int it : {2}) {
  //   double rverse = 0.;
  //   int npt = (it == 2) ? funcs[it].s.m_b.size() : funcs[it].s.m_b.size() -
  //   1; for (int i = 0; i < npt; i++) {
  //     cout << "m_b = " << funcs[it].s.m_b[i] << endl;
  //     rverse += funcs[it].s.m_b[i];
  //   }
  //   if (rverse < 0) return -1;  // return if isn't quaqatic
  // }

  int ls[] = {PHI, RHO};
  for (int it : ls) {
    double invrg = 1. / square11(funcs[it].rng);
    double tm = 0.0;
    for (int i = 0; i < funcs[it].s.m_b.size() - 1; i++)
      tm += (square11(funcs[it].s.m_b[i]) +
             0.5 * funcs[it].s.m_b[i] * funcs[it].s.m_b[i + 1]);
    tm += square11(funcs[it].s.m_b.back());
    error["punish"] += 1e-4 * tm * invrg;
  }

  for (Config& cnf : configs) {
    forceEAM(cnf);
    for (pfAtom& atm : cnf.atoms) {
      for (int it : {X, Y, Z}) {
        atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] - atm.frc[it];
        error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }
    }
    error["frc"] += square11(cnf.fitengy - cnf.engy);
    omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
    ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
  }
  error["shift"] += square11(omaxrho - funcs[EMF].xx.back());
  error["shift"] += square11(ominrho - funcs[EMF].xx.front());
  error["frc"] *= 1e2;
  error["shift"] *= dparams["pshift"] /
                    square11(funcs[EMF].xx.back() - funcs[EMF].xx.front());
  return (error["frc"] * (1 + error["punish"]));  // + error["shift"]);
}

double pfHome::forceEAM(const vector<double>& vv) { return 0.0; }

void pfHome::forceEAM(Config& cnf) {
  cnf.phiengy = cnf.emfengy = 0.0;
  cnf.rhomx = -1e4, cnf.rhomi = 1e4;

  for (pfAtom& atm : cnf.atoms) { /* reset values */
    atm.crho = 0.0;
    for (int it : {X, Y, Z}) atm.phifrc[it] = atm.rhofrc[it] = 0.0;
  }  // ii

  for (pfAtom& atm : cnf.atoms) { /* atoms pairs densities */
    double e0 = funcs[EMF].s(0.0);
    for (Neigh& ngb : atm.neighs) {
      funcs[PHI].s.deriv(ngb.slots[PHI], ngb.shifts[PHI], ngb.phi, ngb.phig);
      cnf.phiengy += ngb.phi;

      double tmp[3];
      for (int it : {X, Y, Z}) {
        tmp[it] = ngb.dist2r[it] * ngb.phig;
        atm.phifrc[it] += tmp[it];
        cnf.atoms[ngb.aid].phifrc[it] -= tmp[it];
      }

      funcs[RHO].s.deriv(ngb.slots[RHO], ngb.shifts[RHO], ngb.rho, ngb.rhog);

      atm.crho += ngb.rho;
      cnf.atoms[ngb.aid].crho += ngb.rho;
    }  // nn

    cnf.rhomx = atm.crho > cnf.rhomx ? atm.crho : cnf.rhomx;
    cnf.rhomi = atm.crho < cnf.rhomi ? atm.crho : cnf.rhomi;

    cnf.emfengy += funcs[EMF].s(atm.crho) - e0;
  }  // ii

  for (pfAtom& atm : cnf.atoms) /* eambedding forces */
    for (Neigh& ngb : atm.neighs) {
      double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
      double tmp[3];
      for (int it : {0, 1, 2}) {
        tmp[it] = ngb.dist2r[it] * emf;
        atm.rhofrc[it] += tmp[it];
        cnf.atoms[ngb.aid].rhofrc[it] -= tmp[it];
      }
    }  // nn
  cnf.fitengy = (cnf.phiengy + cnf.emfengy) / cnf.natoms;
}
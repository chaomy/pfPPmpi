/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-14 23:03:46
 */

#include "pfForce.h"

double pfHome::pfForce::forceEAM(const arma::mat& vv, int tg) {
  while (true) {
    broadcast(cmm, tg, PFROOT);
    if (tg == EXT) break;

    int cnt = 0;
    for (int i : {0, 1, 2}) {
      if (hm.optidx[i] == 0) continue;
      Func& ff = hm.funcs[i];
      for (int j : ff.rlxid) ff.yy[j] = vv[cnt++];
    }

    for (int i = 0; i < hm.funcs.size(); i++) {  // broadcast functions
      broadcast(cmm, hm.funcs[i].xx, PFROOT);
      broadcast(cmm, hm.funcs[i].yy, PFROOT);
    }

    for (Func& ff : hm.funcs) ff.s.set_points(ff.xx, ff.yy);

    hm.error["frc"] = 0.0, hm.error["punish"] = 0.0, hm.error["shift"] = 0.0;
    hm.omaxrho = -1e10, hm.ominrho = 1e10;
    double efrc = 0.0;

    int ls[] = {PHI, RHO};
    for (int it : ls) {
      // double invrg = 1. / square11(hm.funcs[it].rng);
      double tm = 0.0;
      for (int i = 0; i < hm.funcs[it].s.m_b.size() - 1; i++)
        tm += (square11(hm.funcs[it].s.m_b[i]) +
               0.5 * hm.funcs[it].s.m_b[i] * hm.funcs[it].s.m_b[i + 1]);
      tm += square11(hm.funcs[it].s.m_b.back());
      hm.error["punish"] += tm;  //  * invrg;
    }
    hm.error["punish"] = exp(hm.error["punish"] / dparams["pweight"]);

    for (int i = locstt; i < locend; i++) {
      Config& cnf = configs[i];
      forceEAM(cnf);
      for (pfAtom& atm : cnf.atoms) {
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] - atm.frc[it];
          efrc += square11(atm.fitfrc[it] * atm.fweigh[it]);
        }
      }
      efrc += square11(cnf.fitengy - cnf.engy);
      hm.omaxrho = cnf.rhomx > hm.omaxrho ? cnf.rhomx : hm.omaxrho;
      hm.ominrho = cnf.rhomi < hm.ominrho ? cnf.rhomi : hm.ominrho;
    }
    // reduce(cmm, epsf, hm.error["shift"], std::plus<double>(), PFROOT);
    reduce(cmm, efrc, hm.error["frc"], std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  return hm.error["frc"];  // hm.error["punish"] + hm.error["shift"];
}

double pfHome::pfForce::forceEAM(const arma::mat& vv) {
  hm.error["frc"] = 0.0, hm.error["punish"] = 0.0, hm.error["shift"] = 0.0;
  hm.omaxrho = -1e10, hm.ominrho = 1e10;
  int cnt = 0;
  for (int i = 0; i < hm.funcs.size(); i++) { /* interpolates */
    Func& ff = hm.funcs[i];
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

  int ls[] = {PHI, RHO};
  for (int it : ls) {
    double invrg = 1. / square11(hm.funcs[it].rng);
    double tm = 0.0;
    for (int i = 0; i < hm.funcs[it].s.m_b.size() - 1; i++)
      tm += (square11(hm.funcs[it].s.m_b[i]) +
             0.5 * hm.funcs[it].s.m_b[i] * hm.funcs[it].s.m_b[i + 1]);
    tm += square11(hm.funcs[it].s.m_b.back());
    hm.error["punish"] += 1e-4 * tm * invrg;
  }

  for (int i = locstt; i < locend; i++) {
    Config& cnf = configs[i];
    forceEAM(cnf);
    for (pfAtom& atm : cnf.atoms) {
      for (int it : {X, Y, Z}) {
        atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] - atm.frc[it];
        hm.error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }
    }
    hm.error["frc"] += square11(cnf.fitengy - cnf.engy);
    hm.omaxrho = cnf.rhomx > hm.omaxrho ? cnf.rhomx : hm.omaxrho;
    hm.ominrho = cnf.rhomi < hm.ominrho ? cnf.rhomi : hm.ominrho;
  }
  hm.error["shift"] += square11(hm.omaxrho - hm.funcs[EMF].xx.back());
  hm.error["shift"] += square11(hm.ominrho - hm.funcs[EMF].xx.front());
  hm.error["frc"] *= 1e2;
  hm.error["shift"] *= dparams["pshift"] / square11(hm.funcs[EMF].xx.back() -
                                                    hm.funcs[EMF].xx.front());

  double tmp = hm.error["frc"] * (1 + hm.error["punish"]);

  cout << " rank " << cmm.rank() << " err  = " << tmp << endl;
  reduce(cmm, tmp, hm.error["sum"], std::plus<double>(), PFROOT);
  cout << " rank " << cmm.rank() << " err  = " << tmp << endl;
  return hm.error["sum"];
}

void pfHome::pfForce::forceEAM(Config& cnf) {
  cnf.phiengy = cnf.emfengy = 0.0;
  cnf.rhomx = -1e4, cnf.rhomi = 1e4;

  for (pfAtom& atm : cnf.atoms) { /* reset values */
    atm.crho = 0.0;
    for (int it : {X, Y, Z}) atm.phifrc[it] = atm.rhofrc[it] = 0.0;
  }  // ii

  for (pfAtom& atm : cnf.atoms) { /* atoms pairs densities */
    // double e0 = hm.funcs[EMF].s(0.0);
    for (Neigh& ngb : atm.neighs) {
      hm.funcs[PHI].s.deriv(ngb.slots[PHI], ngb.shifts[PHI], ngb.phi, ngb.phig);
      cnf.phiengy += ngb.phi;

      double tmp[3];
      for (int it : {X, Y, Z}) {
        tmp[it] = ngb.dist2r[it] * ngb.phig;
        atm.phifrc[it] += tmp[it];
        cnf.atoms[ngb.aid].phifrc[it] -= tmp[it];
      }

      hm.funcs[RHO].s.deriv(ngb.slots[RHO], ngb.shifts[RHO], ngb.rho, ngb.rhog);

      atm.crho += ngb.rho;
      cnf.atoms[ngb.aid].crho += ngb.rho;
    }  // nn

    double embE;
    hm.funcs[EMF].s.deriv(atm.crho, embE, atm.gradF);
    cnf.emfengy += (embE);

    cnf.rhomx = atm.crho > cnf.rhomx ? atm.crho : cnf.rhomx;
    cnf.rhomi = atm.crho < cnf.rhomi ? atm.crho : cnf.rhomi;
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
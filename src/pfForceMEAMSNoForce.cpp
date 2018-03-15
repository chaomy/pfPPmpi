/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-15 16:41:18
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

double pfHome::forceMEAMSNoForce(const arma::mat &vv, int tg) {
  while (true) {
    broadcast(cmm, tg, PFROOT);
    if (tg == EXT) break;

    int cnt = 0;
    for (int i : {0, 1, 2, 3, 4}) {
      if (optidx[i] == 0) continue;
      Func &ff = funcs[i];
      for (int j : ff.rlxid) ff.yy[j] = vv[cnt++];
    }

    for (int i = 0; i < nfuncs; i++) {  // broadcast functions
      broadcast(cmm, funcs[i].xx, PFROOT);
      broadcast(cmm, funcs[i].yy, PFROOT);
    }

    for (Func &ff : funcs) ff.s.set_points(ff.xx, ff.yy);

    double efrc = 0.0, eengy = 0.0;
    error["punish"] = 0.0;
    omaxrho = -1e10, ominrho = 1e10;

    int ww = 1;
    for (int it : smthidx) {  // covarance of third derivative
      vector<double> &vv = funcs[it].s.m_a;
      double mn = 0.0, cov = 0.0;
      for (int i = ww + 1; i < vv.size() - ww; i++) {
        for (int it = -ww; it <= ww; it++) mn += vv[i + it];
        mn /= (2 * ww + 1);
        for (int it = -ww; it <= ww; it++) cov += square11(vv[i + it] - mn);
      }
      error["punish"] += cov / (2 * ww + 1);
      // error["punish"] += square11(mn);
    }
    error["punish"] *= dparams["pweight"];

    for (int i : locls) {
      Config &cnf = configs[i];
      forceMEAMS(cnf);
      eengy += cnf.weigh * square11(cnf.fitengy - cnf.engy);
      omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
      ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
    }

    eengy *= dparams["eweight"];
    reduce(cmm, eengy, error["engy"], std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  return error["engy"] + error["punish"];
}

void pfHome::forceMEAMSNoForce(Config &cnf) {  // It's benchmark one
  cnf.phiengy = cnf.emfengy = 0.0;
  cnf.rhomi = 1e10, cnf.rhomx = -1e10;
  for (pfAtom &atm : cnf.atoms) atm.crho = 0.0;

  // ii
  double e0 = funcs[EMF].s(0.0);
  for (pfAtom &atm : cnf.atoms) { /* loop over atoms pairs, densities */
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];

      // pair (phi)
      funcs[PHI].s.deriv(ngbj.slots[PHI], ngbj.shifts[PHI], ngbj.phi,
                         ngbj.phig);
      cnf.phiengy += ngbj.phi;
      for (int it : {X, Y, Z}) atm.phifrc[it] += ngbj.dist2r[it] * ngbj.phig;

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
  }                             // atm
  for (pfAtom &atm : cnf.atoms) /* eambedding forces */
    for (Neigh &ngb : atm.neighsFull) {
      double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
      for (int it : {X, Y, Z}) atm.rhofrc[it] += ngb.dist2r[it] * emf;
    }  // nn
  cnf.fitengy = (cnf.phiengy / 2. + cnf.emfengy) / cnf.natoms;
}

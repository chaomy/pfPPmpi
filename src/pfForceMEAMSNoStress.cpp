/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-19 14:01:38
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

double pfHome::forceMEAMS(const arma::mat &vv, int tg) {
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
    error["frc"] = 0.0, error["punish"] = 0.0;
    double omax = -1e10, omin = 1e10;

    int ww = 1;

    // to inference covarance of third derivative
    // for (int it : smthidx) {
    //   vector<double> &vv = funcs[it].s.m_a;
    //   double mn = 0.0, cov = 0.0;
    //   for (int i = ww + 1; i < vv.size() - ww; i++) {
    //     for (int it = -ww; it <= ww; it++) mn += vv[i + it];
    //     mn /= (2 * ww + 1);
    //     for (int it = -ww; it <= ww; it++) cov += square11(vv[i + it] - mn);
    //   }
    //   error["punish"] += cov / (2 * ww + 1);
    //   // error["punish"] += square11(mn);
    // }
    // error["punish"] *= dparams["pweight"];
    // to decrease third derivative
    // for (int it : smthidx) {
    //   vector<double> &vv = funcs[it].s.m_a;
    //   for (auto ee : vv) error["punish"] += square11(ee);
    // }

    double rs = 0;
    double Mf = dparams["fwidth"], Me = dparams["ewidth"];
    double Bf = dparams["lrange"] * Mf, Be = dparams["lrange"] * Me;
    double OutE = Me * (2 * Be - Me);
    double OutF = Mf * (2 * Bf - Mf);

    for (int i : locls) {
      Config &cnf = configs[i];
      forceMEAMS(cnf);
      for (pfAtom &atm : cnf.atoms) {
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] =
              atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it] - atm.frc[it];
          rs = fabs(atm.fitfrc[it] * atm.fweigh[it]);
          if (rs < Bf)
            efrc += cnf.weigh * (rs < Mf ? square11(rs) : Mf * (2 * rs - Mf));
          else
            efrc += cnf.weigh * OutF;
        }
      }
      rs = fabs(cnf.fitengy - cnf.engy);
      if (rs < Be)
        eengy += cnf.weigh * (rs < Me ? square11(rs) : Me * (2 * rs - Me));
      else
        eengy += cnf.weigh * OutE;
      omax = cnf.rhomx > omax ? cnf.rhomx : omax;
      omin = cnf.rhomi < omin ? cnf.rhomi : omin;
    }
    eengy *= dparams["eweight"];
    reduce(cmm, omin, ominrho, mpi::minimum<double>(), PFROOT);
    reduce(cmm, omax, omaxrho, mpi::maximum<double>(), PFROOT);
    reduce(cmm, eengy, error["engy"], std::plus<double>(), PFROOT);
    reduce(cmm, efrc, error["frc"], std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  error["engy"] *= 100;
  error["frc"] *= 100;
  return error["frc"] + error["engy"] + error["punish"];
}

void pfHome::forceMEAMS(Config &cnf) {  // It's benchmark one
  cnf.phiengy = cnf.emfengy = 0.0;
  cnf.rhomi = 1e10, cnf.rhomx = -1e10;
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
      for (int it : {X, Y, Z}) atm.rhofrc[it] += ngb.dist2r[it] * emf;
      // }  // cutoff(rho)
    }  // nn
  cnf.fitengy = (cnf.phiengy / 2. + cnf.emfengy) / cnf.natoms;
}

// close fiting physical
// error["phy"] = 0.0;
// if (iparams["runlmp"]) {
//   (this->*write[sparams["ptype"]])();
//   lmpdrv->calLatticeBCC();
//   lmpdrv->calLatticeFCC();
//   lmpdrv->calLatticeHCP();
//   lmpdrv->calSurfaceUrlx();
//   lmpdrv->calGSFUrlx();

//   remove("no");
//   remove("log.lammps");
//   remove("restart.equil");
//   lmpdrv->exprs["bcc2hcp"] = lmpdrv->exprs["ehcp"] - lmpdrv->exprs["ebcc"];
//   lmpdrv->exprs["bcc2fcc"] = lmpdrv->exprs["efcc"] - lmpdrv->exprs["ebcc"];

//   vector<string> aa(
//       {"lat", "bcc2fcc", "bcc2hcp", "suf110", "suf100", "suf111"});
//   vector<double> ww({2e5, 5e3, 5e3, 1e3, 1e3, 1e3});

//   for (int i = 0; i < aa.size(); i++) {
//     string ee(aa[i]);
//     error["phy"] +=
//         (lmpdrv->error[ee] =
//              ww[i] * square11((lmpdrv->exprs[ee] - lmpdrv->targs[ee]) /
//                               lmpdrv->targs[ee]));
//   }

//   error["gsf"] = 0.0;
//   for (int i : lmpdrv->gsfpnts)
//     error["gsf"] +=
//         500 * (lmpdrv->lgsf["111e110"][i] + lmpdrv->lgsf["111e211"][i]);
//   error["phy"] += error["gsf"];
//   error["phy"] *= phyweigh;
// }

// int ls[] = {PHI, RHO, MEAMF};
// for (int it : ls) {
//   double invrg = 1. / square11(funcs[it].rng);
//   double tm = 0.0;
//   for (int i = 0; i < funcs[it].s.m_b.size() - 1; i++)
//     tm += (square11(funcs[it].s.m_b[i]) +
//            0.5 * funcs[it].s.m_b[i] * funcs[it].s.m_b[i + 1]);
//   tm += square11(funcs[it].s.m_b.back());
//   epsh += tm * invrg;
// }
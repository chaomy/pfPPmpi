/*
 * @Author: chaomy
 * @Date:   2017-12-13 09:53:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-19 15:21:04
 */

#include "pfHome.h"

using std::setw;

void pfHome::resample() {
  for (auto kk : {0, 1, 2, 3}) {
    int npts = funcs[kk].npts;
    // double ri = ricut - 0.01;
    double ri = funcs[PHI].xx[0];
    double ro = 0.0;
    if (kk == PHI)
      ro = 5.25;
    else if (kk == RHO || kk == MEAMF)
      ro = 5.25;
    else if (kk == EMF) {
      ri = ominrho - 40.0;
      ro = omaxrho + 20.0;
      double dl = (funcs[EMF].xx.back() - funcs[EMF].xx.front()) / 3.;
      for (int i = 0; i <= 3; i++) {
        double xx = dl * i + funcs[EMF].xx.front();
        double phinew, phinewdriv;
        funcs[EMF].s.deriv(xx, phinew, phinewdriv);
        cout << std::setprecision(16) << xx << " " << phinew << " "
             << phinewdriv << endl;
      }
    } else if (kk == MEAMG) {
      ri = -1.0;
      ro = 1.0;
    }
    double delt = (ro - ri) / (npts - 1);
    for (int i = 0; i < npts; i++) {
      funcs[kk].xx[i] = delt * i + ri;
      funcs[kk].s.deriv(funcs[kk].xx[i], funcs[kk].yy[i], funcs[kk].g1[i]);
    }
  }
  for (Func& ff : funcs) ff.s.set_points(ff.xx, ff.yy);
  (this->*write[sparams["ptype"]])();
}

void pfHome::analyLoss() {
  double rs = 0;
  double Mf = dparams["fbndq"], Me = dparams["ebndq"];
  double Bf = dparams["fbndl"], Be = dparams["ebndl"];
  double OutE = Me * (2 * Be - Me);
  double OutF = Mf * (2 * Bf - Mf);

  int qEnb = 0, lEnb = 0, cEnb = 0;
  int qFnb = 0, lFnb = 0, cFnb = 0;
  double qEer = 0, lEer = 0, cEer = 0;
  double qFer = 0, lFer = 0, cFer = 0;

  for (Config& cnf : configs) {
    for (pfAtom& atm : cnf.atoms) {
      for (int it : {X, Y, Z}) {
        rs = fabs(atm.fitfrc[it] * atm.fweigh[it]);
        if (rs < Mf) {
          qFnb += 1;
          qFer += square11(rs);
        } else if (rs < Bf) {
          lFnb += 1;
          lFer += Mf * (2 * rs - Mf);
        } else {
          cFnb += 1;
          cFer += OutF;
        }
      }
    }
    rs = fabs(cnf.fitengy - cnf.engy);
    if (rs < Me) {
      qEnb += 1;
      qEer += square11(rs);
    } else if (rs < Be) {
      lEnb += 1;
      lEer += Me * (2 * rs - Me);
    } else {
      cEnb += 1;
      cEer += OutE;
    }
  }
  cout << setw(10) << "Enum " << qEnb << " " << lEnb << " " << cEnb << endl;
  cout << setw(10) << "Eerr " << qEer << " " << lEer << " " << cEer << endl;
  cout << setw(10) << "Fnum " << qFnb << " " << lFnb << " " << cFnb << endl;
  cout << setw(10) << "Ferr " << qFer << " " << lFer << " " << cFer << endl;
}

void pfHome::writeRadDist() {
  ofstream of("radius.txt", std::ofstream::out);
  for (auto& cc : configs)
    for (auto& atm : cc.atoms)
      for (auto& ngb : atm.neighsFull) {
        of << ngb.r << " " << funcs[PHI].s(ngb.r) << " " << funcs[RHO].s(ngb.r)
           << " " << funcs[MEAMF].s(ngb.r) << endl;
      }
  of.close();
}

void pfHome::writeAngDist() {
  ofstream of("angle.txt", std::ofstream::out);
  for (auto& cc : configs) {
    for (auto& atm : cc.atoms) {
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        for (int kk = 0; kk < jj; kk++)
          of << atm.angMat[jj][kk].gcos << " " << atm.angMat[jj][kk].gval << " "
             << atm.angMat[jj][kk].ggrad << endl;
      }
    }
  }
  of.close();
}

void pfHome::writeFrcDist() {  // check the behaviors of force fitting
  double M = dparams["fwidth"];
  ofstream of("force.txt", std::ofstream::out);
  for (int i : locls) {
    for (pfAtom& atm : configs[i].atoms)
      for (int it : {0, 1, 2}) {
        double rs = fabs(atm.fitfrc[it] * atm.fweigh[it]);
        double err = rs < M ? square11(rs) : M * (2 * rs - M);
        of << std::setprecision(6) << atm.frc[it] << " "
           << atm.fitfrc[it] + atm.frc[it] << " " << atm.phifrc[it] << " "
           << atm.rhofrc[it] << " " << atm.trifrc[it] << " " << rs << " "
           << square11(rs) << " " << err << endl;
      }
  }
  of.close();
}

void pfHome::writeEngDist() {  // check energy distributions
  double M = dparams["ewidth"];
  ofstream of("engy.txt", std::ofstream::out);
  for (auto& cnf : configs) {
    double rs = fabs(cnf.fitengy - cnf.engy);
    double err = rs < M ? square11(cnf.fitengy - cnf.engy) : M * (2 * rs - M);
    of << std::setprecision(6) << cnf.fitengy << " " << cnf.engy << " " << rs
       << " " << square11(cnf.fitengy - cnf.engy) << " " << err << endl;
  }
  of.close();
}

void pfHome::loopBwth() {
  for (int i = 1; i < 20; i++) {
    double tm = 0.2 + 0.2 * i;
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        for (int it : {0, 1, 2})
          atm.fweigh[it] = exp(-tm * atm.frc[it] * atm.frc[it]);

    vector<double> mtm(3, 0);
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        for (int ii : {0, 1, 2}) mtm[ii] += fabs(atm.frc[ii]) * atm.fweigh[ii];

    for (int ii : {0, 1, 2}) mtm[ii] /= ftn;

    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        for (int ii : {0, 1, 2}) atm.fweigh[ii] *= (mfrc[ii] / mtm[ii]);

    string fm("fc" + to_string(i) + ".txt");
    FILE* fid = fopen(fm.c_str(), "w");
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        fprintf(fid, "%f %f %f %f %f %f\n", atm.frc[X], atm.frc[Y], atm.frc[Z],
                atm.frc[X] * atm.fweigh[X], atm.frc[Y] * atm.fweigh[Y],
                atm.frc[Z] * atm.fweigh[Z]);
    fclose(fid);
  }  // i = 1 to 9
}

void pfHome::forceDis() {
  vector<double> mtm(3, 0);
  for (Config& tmpc : configs)
    for (pfAtom& atm : tmpc.atoms)
      for (int ii : {0, 1, 2}) mtm[ii] += fabs(atm.frc[ii]) * atm.fweigh[ii];
  fsm = square11(mfrc[X]) + square11(mfrc[Y]) + square11(mfrc[Z]);
}

void pfHome::deleteAtoms() { /* delete some atoms (Better Not)*/
  for (Config& tmpc : configs) {
    vector<pfAtom> ntm;
    for (pfAtom& atm : tmpc.atoms)
      if (atm.nneighs == N3 && atm.nneighsFull == N3) ntm.push_back(atm);
    tmpc.atoms = ntm;
    tmpc.natoms = tmpc.atoms.size();
  }
}

void pfHome::cutoffNeighs() { /* for chossing a cutoff*/
  FILE* fid = fopen("cutoff.txt", "w");
  rocut = 4.60;

  /* loop over cutoff */
  while (rocut < 6.40) {
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms) atm.neighsFull.clear();

    initNeighsFull();
    fprintf(fid, "%0.3f ", rocut);
    for (Config& tmpc : configs) {
      int aven = 0;
      for (pfAtom& atm : tmpc.atoms) aven += atm.nneighsFull;
      fprintf(fid, "%03.1f ", float(aven) / tmpc.natoms);
    }

    fprintf(fid, "\n");
    rocut += 0.01;
  }

  fclose(fid);
}

// void pfHome::resample() {
//   int npts = 12;
//   double rout = 6.25;
//   ricut -= 0.01;
//   double delt = (rout - ricut) / npts;
//   cout << "PHI" << endl;

//   for (int i = 0; i <= npts; i++) {
//     double xx = delt * i + ricut;
//     double phinew, phinewdriv;
//     funcs[PHI].s.deriv(xx, phinew, phinewdriv);
//     cout << std::setprecision(16) << xx << " " << phinew << " " << phinewdriv
//          << endl;
//   }

//   npts = 11;
//   rout = 5.50;
//   delt = (rout - ricut) / npts;
//   cout << "RHO" << endl;
//   for (int i = 0; i <= npts; i++) {
//     double xx = delt * i + ricut;
//     double phinew, phinewdriv;
//     funcs[RHO].s.deriv(xx, phinew, phinewdriv);
//     cout << std::setprecision(16) << xx << " " << phinew << " " << phinewdriv
//          << endl;
//   }

//   npts = 11;
//   rout = 5.50;
//   delt = (rout - ricut) / npts;
//   cout << "MEAMF" << endl;
//   for (int i = 0; i <= npts; i++) {
//     double xx = delt * i + ricut;
//     double phinew, phinewdriv;
//     funcs[MEAMF].s.deriv(xx, phinew, phinewdriv);
//     cout << std::setprecision(16) << xx << " " << phinew << " " << phinewdriv
//          << endl;
//   }

//   npts = 3;
//   double rhmn = -100.0;
//   double rhmx = -60.0;
//   double dl = (rhmx - rhmn) / (npts - 1);

//   cout << "RHO" << endl;
//   for (int i = 0; i <= npts; i++) {
//     double xx = dl * i + rhmn;
//     double phinew, phinewdriv;
//     funcs[EMF].s.deriv(xx, phinew, phinewdriv);
//     cout << std::setprecision(16) << xx << " " << phinew << " " << phinewdriv
//          << endl;
//   }

//   dl = ((rhmx = 1. + 1e-12) - (rhmn = -1.0 - 1e-12)) / ((npts = 6) - 1);
//   cout << "MEAMG" << endl;
//   for (int i = 0; i < npts; i++) {
//     double xx = dl * i + rhmn;
//     double phinew, phinewdriv;
//     funcs[MEAMG].s.deriv(xx, phinew, phinewdriv);
//     cout << std::setprecision(16) << xx << " " << phinew << " " << phinewdriv
//          << endl;
//   }
// }
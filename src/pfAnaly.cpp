/*
 * @Author: chaomy
 * @Date:   2017-12-13 09:53:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-16 15:56:29
 */

#include "pfHome.h"

using std::cout;
using std::endl;
using std::to_string;
using std::vector;

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
      ri = ominrho - 15.0;
      ro = omaxrho + 15.0;
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
  for (Config& tmpc : configs) {
    for (pfAtom& atm : tmpc.atoms)
      for (int ii : {0, 1, 2}) mtm[ii] += fabs(atm.frc[ii]) * atm.fweigh[ii];
  }
  // for (int ii : {0, 1, 2}) mtm[ii] /= ftn;
  // for (Config& tmpc : configs) {
  //   for (pfAtom& atm : tmpc.atoms)
  //     for (int ii : {0, 1, 2}) atm.fweigh[ii] *= (mfrc[ii] / mtm[ii]);
  // }
  fsm = square11(mfrc[X]) + square11(mfrc[Y]) + square11(mfrc[Z]);
}

void pfHome::deleteAtoms() { /* delete some atoms */
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
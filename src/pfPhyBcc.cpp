/*
 * @Author: chaomy
 * @Date:   2017-12-16 16:00:09
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-26 16:09:37
 */

#include "pfConf.h"
#include "pfForce.h"
#include "pfPhy.h"

void pfHome::pfPhy::calLat(string kk, int npts, pfForce& fcdrv, pfConf& cdrv) {
  double la = (kk == "bcc") ? 3.30 : 4.20;
  double del = 1.5;
  double dl = del / npts;
  double lo = la - del;

  ofstream ostr((kk + ".txt"), std::ofstream::out);
  for (int i = 0; i < 2 * npts + 1; i++) {
    double latt = lo + i * dl;
    // Config cc(buildbccPrim(latt));
    Config cc = (cdrv.*cdrv.build[kk])(latt);
    (fcdrv.*calfrc[sparams["ptype"]])(cc);
    ostr << latt << " " << cc.atoms[0].crho << " " << cc.fitengy << " "
         << 0.5 * cc.phiengy / cc.natoms << " " << cc.emfengy / cc.natoms
         << endl;
  }
  ostr.close();
}

void pfHome::pfPhy::calLat(string kk, pfForce& fcdrv, pfConf& cdrv) {
  double la = 3.30, me = 1e3;
  unordered_map<string, vector<double>>& mpvc = cdrv.mpvc;
  for (double dl : {3e-2, 9e-3, 3e-3, 9e-4, 3e-4}) {
    if (!kk.compare("bcc")) cdrv.buildbcc(kk, la, dl);
    for (int ii = 0; ii < mpvc[kk].size(); ii++) {
      Config& cc = cdrv.mpcf[kk][ii];
      (fcdrv.*calfrc[sparams["ptype"]])(cc);
      if (cc.fitengy < me) {
        me = cc.fitengy;
        la = mpvc[kk][ii];
      }
    }
  }
  if (!kk.compare("bcc")) cdrv.ubcc = cdrv.buildbccPrim(la);
  exprs["lat"] = la;
  perr["lat"] = pwgh["lat"] * square11(la - targs["lat"]);
  (fcdrv.*calfrc[sparams["ptype"]])(cdrv.ubcc);
}

void pfHome::pfConf::buildbcc(const string& kk, const double& gs,
                              const double& dl) {
  double l = gs - 3 * dl;
  mpvc[kk].clear();
  mpcf[kk].clear();
  for (int i = 0; i < 7; i++) {
    mpvc[kk].push_back(l + i * dl);
    mpcf[kk].push_back(buildbccPrim(mpvc[kk].back()));
  }
}

Config pfHome::pfConf::buildbccPrim(const double& lat) {
  Config cc;
  double ll = 0.5 * lat;
  cc.bvx[X] = -ll, cc.bvx[Y] = ll, cc.bvx[Z] = ll;
  cc.bvy[X] = ll, cc.bvy[Y] = -ll, cc.bvy[Z] = ll;
  cc.bvz[X] = ll, cc.bvz[Y] = ll, cc.bvz[Z] = -ll;

  pfAtom atm(0);
  atm.prl[X] = atm.prl[Y] = atm.prl[Z] = 0.0;
  atm.pst[X] = atm.pst[Y] = atm.pst[Z] = 0.0;
  cc.atoms.push_back(atm);

  cc.natoms = cc.atoms.size();
  initBox(cc);
  if (!sparams["ptype"].compare("MEAMS")) {
    initNeighsFull(cc);
    initAngles(cc);
  } else if (!sparams["ptype"].compare("EAMS"))
    initNeighs(cc);
  return cc;
}

Config pfHome::pfConf::buildbccConv(const double& lat) {
  Config cc;  // use 2 x 2 x 2
  cc.bvx[X] = lat, cc.bvx[Y] = 0.0, cc.bvx[Z] = 0.0;
  cc.bvy[X] = 0.0, cc.bvy[Y] = lat, cc.bvy[Z] = 0.0;
  cc.bvz[X] = 0.0, cc.bvz[Y] = 0.0, cc.bvz[Z] = lat;
  int cnt = 0;

  // add atoms
  pfAtom atm1(cnt++);
  atm1.prl[X] = atm1.prl[Y] = atm1.prl[Z] = 0.0;
  for (int i : {X, Y, Z})
    atm1.pst[i] = atm1.prl[X] * cc.bvx[i] + atm1.prl[Y] * cc.bvy[i] +
                  atm1.prl[Z] * cc.bvz[i];
  cc.atoms.push_back(atm1);

  pfAtom atm2(cnt++);
  atm2.prl[X] = atm2.prl[Y] = atm2.prl[Z] = 0.5;
  for (int i : {X, Y, Z})
    atm2.pst[i] = atm2.prl[X] * cc.bvx[i] + atm2.prl[Y] * cc.bvy[i] +
                  atm2.prl[Z] * cc.bvz[i];
  cc.atoms.push_back(atm2);

  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}
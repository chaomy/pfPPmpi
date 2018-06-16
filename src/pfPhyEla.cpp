/*
 * @Author: chaomy
 * @Date:   2017-12-17 11:13:45
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-16 00:32:11
 */

#include "pfConf.h"
#include "pfForce.h"
#include "pfPhy.h"

#define C11A (1. / 9.)
#define C11B (1. / 3.)
#define C12A (1. / 9.)
#define C12B (1. / 6.)

void pfHome::pfPhy::calElas(int npts, pfForce& fcdrv, pfConf& cfdrv) {
  Config c1(cfdrv.buildbccPrim(exprs["lat"]));
  double delta = 0.01;
  double lo = -delta * npts;

  ofstream ostr("c12.txt", std::ofstream::out);
  for (int i = 0; i < 2 * npts + 1; i++) {
    double dl = lo + delta * i;

    Config cc = cfdrv.addotho(c1, dl);
    (fcdrv.*calfrc[sparams["ptype"]])(cc);
    ostr << dl << " " << cc.atoms[0].crho << " " << cc.fitengy << " "
         << 0.5 * cc.phiengy / cc.natoms << " " << cc.emfengy / cc.natoms
         << endl;
  }
  ostr.close();
}

void pfHome::pfPhy::calElas(pfForce& fcdrv, pfConf& cfdrv) {
  vector<double> dv({-3e-3, -9e-4, -3e-4, 3e-4, 9e-4, 3e-3});
  unordered_map<string, vector<Config>>& mpcf = cfdrv.mpcf;
  Config& ubcc = cfdrv.ubcc;

  for (string kk : {"c11", "c12", "c44"}) cfdrv.mpcf[kk].clear();
  for (double dl : dv) {
    mpcf["c11"].push_back(cfdrv.addvolm(ubcc, dl));
    mpcf["c12"].push_back(cfdrv.addotho(ubcc, dl));
    mpcf["c44"].push_back(cfdrv.addmono(ubcc, dl));
  }

  for (auto& ee : mpcf["c11"]) (fcdrv.*calfrc[sparams["ptype"]])(ee);
  for (auto& ee : mpcf["c12"]) (fcdrv.*calfrc[sparams["ptype"]])(ee);
  for (auto& ee : mpcf["c44"]) (fcdrv.*calfrc[sparams["ptype"]])(ee);

  vector<double> tm(3, 0);
  // double vol = 0.5 * exprs["lat"] * exprs["lat"] * exprs["lat"];
  double vol = ubcc.vol;

  for (int i = 0; i < dv.size(); i++) {
    tm[0] += 2. * (mpcf["c11"][i].fitengy - ubcc.fitengy) / (dv[i] * dv[i]);
    tm[1] += 2. * (mpcf["c12"][i].fitengy - ubcc.fitengy) / (dv[i] * dv[i]);
    tm[2] += 2. * (mpcf["c44"][i].fitengy - ubcc.fitengy) / (dv[i] * dv[i]);
  }

  for (int i : {0, 1, 2}) tm[i] /= dv.size();
  for (int i : {0, 1, 2}) tm[i] *= (EVA3_GPA / vol);

  exprs["c11"] = C11A * tm[0] + C11B * tm[1];
  exprs["c12"] = C12A * tm[0] - C12B * tm[1];
  exprs["c44"] = 0.25 * tm[2];

  error["ela"] = 0.0;
  for (auto kk : {"c11", "c12", "c44"})
    error["ela"] += square11(targs[kk] - exprs[kk]) * weigh[kk];
}
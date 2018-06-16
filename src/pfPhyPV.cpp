/*
 * @Author: chaomy
 * @Date:   2017-12-19 08:58:15
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 22:18:30
 */

#include "pfConf.h"
#include "pfPhy.h"
#include "pfForce.h"

using std::pow;
using std::to_string;

void pfHome::pfPhy::calPV(pfForce& fcdrv, pfConf& cfdrv) {
  string ptg(sparams["elem"] + "p");
  string vtg(sparams["elem"] + "v");
  string etg(sparams["elem"] + "e");
  unordered_map<string, vector<Config>>& mpcf = cfdrv.mpcf;

  const vector<double>& pp = mele.pvm[ptg];
  const vector<double>& vv = mele.pvm[vtg];

  mpcf["pv"].clear();
  for (double dl : vv)
    mpcf["pv"].push_back(cfdrv.addvolm(cfdrv.ubcc, pow(dl, 1. / 3.) - 1.0));

  error["pv"] = 0.0;
  for (int i = 0; i < mpcf["pv"].size(); i++) {
    auto& ee = mpcf["pv"][i];
    fcdrv.stressMEAM(ee);
    for (int i = 0; i < 6; i++) ee.strs[i] *= EVA3_GPA;
    error["pv"] += square11(pp[i] + ee.strs[XX]);
  }
}
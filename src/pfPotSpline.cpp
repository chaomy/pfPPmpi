/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-08 22:43:03
 */

#include "pfHome.h"

/* set boundary conditions the lo, hi and delta of the node values */

void pfHome::setSplineVariables() {  // func -> ini
  nvars = 0;
  ini.clear(), lob.clear(), hib.clear(), deb.clear();
  for (Func& ff : funcs) ff.step = ff.xx[1] - ff.xx[0];
  for (int i : {0, 1, 2, 3, 4}) {
    if (optidx[i] == 0) continue;
    Func& ff = funcs[i];
    startps.push_back(nvars);
    double mn = 1e30, mx = -1e30;
    for (int j : ff.rlxid) {
      ini.push_back(ff.yy[j]);
      mn = std::min(ff.yy[j], mn);
      mx = std::max(ff.yy[j], mx);
      nvars++;
    }
    double del = 2 * (mx - mn);  // be careful on the delta
    for (int j : ff.rlxid) {
      lob.push_back(ff.yy[j] - del);
      hib.push_back(ff.yy[j] + del);
      deb.push_back(2 * del);
    }
    endps.push_back(nvars);
  }  // i
}
/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-14 03:43:46
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

// void pfHome::setSplineVariablesNear() {  // func -> ini
//   nvars = 0;
//   ini.clear();
//   lob.clear();
//   hib.clear();
//   deb.clear();
//   for (Func& ff : funcs) ff.step = ff.xx[1] - ff.xx[0];
//   for (int i : {0, 1, 2, 3, 4}) {
//     if (optidx[i] == 0) continue;
//     Func& ff = funcs[i];
//     startps.push_back(nvars);
//     for (int j : ff.rlxid) {
//       ini.push_back(ff.yy[j]);
//       double vari =
//           (fabs(ff.yy[j]) >= 1e-8) ? dparams["ivari"] * fabs(ff.yy[j]) :
//           0.001;
//       lob.push_back(ff.yy[j] - vari);
//       hib.push_back(ff.yy[j] + vari);
//       deb.push_back(2 * vari);
//       nvars++;
//     }
//     endps.push_back(nvars);
//   }  // i
// }

void pfHome::setSplineVariables() {  // func -> ini
  nvars = 0;
  ini.clear();
  lob.clear();
  hib.clear();
  deb.clear();
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
    double del = (mx - mn);  // be careful on the delta
    for (int j : ff.rlxid) {
      lob.push_back(ff.yy[j] - del);
      hib.push_back(ff.yy[j] + del);
      deb.push_back(2 * del);
    }
    endps.push_back(nvars);
  }  // i
}
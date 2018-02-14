/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-13 12:50:31
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

void pfHome::setSplineVariables() {
  // func -> ini
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
    // for (int j = 0; j < ff.npts; j++) {
    for (int j : ff.rlxid) {
      cout << "relx id = " << j << endl;
      ini.push_back(ff.yy[j]);
      double vari =
          (fabs(ff.yy[j]) >= 1e-8) ? dparams["ivari"] * fabs(ff.yy[j]) : 0.001;
      lob.push_back(ff.yy[j] - vari);
      hib.push_back(ff.yy[j] + vari);
      deb.push_back(2 * vari);
      nvars++;
    }
    endps.push_back(nvars);
  }  // i
}
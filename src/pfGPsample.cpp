/*
 * @Author: chaomy
 * @Date:   2018-02-12 21:51:52
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 15:28:00
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
using std::cout;
using std::endl;

/* This function has not been finished yet */

void pfHome::GPsample() {
  ofstream of1("err.txt", std::ofstream::out);
  ofstream of2("par.txt", std::ofstream::out);
  double overallphy = 1e30;

  for (int it = 0; it < iparams["maxstep"]; it++) {
    vector<double> vv =
        decodestdv(5.0 + dparams["istep"] *
                             (arma::mat(nvars, 1, arma::fill::randu) - 0.5));
    cout << "start sample it : " << it << endl;

    of2 << it;
    for (int i = 0; i < nvars; i++) of2 << " " << vv[i];
    of2 << endl;

    int cnt = 0;
    for (int i : {0, 1, 2, 3, 4}) {
      if (optidx[i] == 0) continue;
      Func& ff = funcs[i];
      for (int j : ff.rlxid) ff.yy[j] = vv[cnt++];
    }

    for (Func& ff : funcs) ff.s.set_points(ff.xx, ff.yy);

    (this->*write[sparams["ptype"]])();
    lmpCheck(it, of1);
    if (error["phy"] < overallphy) {
      std::rename(sparams["lmpfile"].c_str(), "meam.lib.best");
      overallphy = error["phy"];
    }
  }

  of1.close();
  of2.close();
}
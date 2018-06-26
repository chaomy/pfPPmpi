/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:31:59
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-20 14:51:35
 */

#include "pfForce.h"
#include "pfConf.h"
#include "pfIO.h"
#include "pfLmpDrv.h"
#include "pfOptimizer.h"
#include "pfPhy.h"

void pfHome::increAnneal(pfForce &fcdrv, pfConf &cfdrv, pfIO &io) {
  scnt = 0;
  int jobl[] = {PHI, PHI, RHO, MEAMF, MEAMG, PHI, RHO, PHI, RHO, MEAMF, MEAMG};
  for (int i = 0; i < 10; i++) {
    simAnneal(fcdrv, io);
    upgrade(jobl[i], fcdrv, cfdrv);
    iparams["kmax"]++;
  }
}

void pfHome::run(int argc, char *argv[], pfForce &fcdrv, pfConf &cfdrv,
                 pfPhy &phdrv, pfIO &io) {
  if (!sparams["opt"].compare("gp"))
    GPsample(io);
  else if (!sparams["opt"].compare("anneal"))
    simAnneal(fcdrv, io);
  else if (!sparams["opt"].compare("evo"))
    diffEvo(fcdrv, io);
  else if (!sparams["opt"].compare("cmaes"))
    loopcmaes(fcdrv, io);
  else if (!sparams["opt"].compare("cnt"))
    cntcmaes(fcdrv, io);
  else if (!sparams["opt"].compare("incre"))
    increAnneal(fcdrv, cfdrv, io);
  else if (!sparams["opt"].compare("nlopt"))
    nloptGlobal(io);
  else if (!sparams["opt"].compare("shift"))
    doShift(fcdrv, io);
  else if (!sparams["opt"].compare("buildD03"))
    cfdrv.buildD03("d03", 7.400, 0.005);
  else if (!sparams["opt"].compare("make")) {
    calErr(fcdrv, io);
    resample(io);
  } else if (!sparams["opt"].compare("anlz")) {
    calErr(fcdrv, io);
    phdrv.calLat("bcc", 30, fcdrv, cfdrv);
    phdrv.calLat("fcc", 30, fcdrv, cfdrv);
    exprs["lat"] = exprs["lat"];
    phdrv.calElas(25, fcdrv, cfdrv);

    exprs["bcc2hcp"] = exprs["ehcp"] - exprs["ebcc"];
    exprs["bcc2fcc"] = exprs["efcc"] - exprs["ebcc"];

    vector<string> aa(
        {"lat", "bcc2fcc", "bcc2hcp", "suf110", "suf100", "suf111"});
    for (auto ee : aa) cout << ee << " " << exprs[ee] << endl;
  }
}

void pfHome::calErr(pfForce &fcdrv, pfIO &io) {  // make potential
  arma::mat mm(nvars, 1, arma::fill::randu);
  for (int i = 0; i < nvars; i++) mm[i] = ini[i];
  (fcdrv.*calobj[sparams["ptype"]])(mm, 1);
  if (cmm.rank() == PFROOT) {
    double err = (fcdrv.*calobj[sparams["ptype"]])(mm, 1);
    (fcdrv.*calobj[sparams["ptype"]])(mm, EXT);
    (io.*write[sparams["ptype"]])();

    ofstream of("err.txt", std::ofstream::out);
    of << std::setprecision(4) << err << " " << error["frc"] << " "
       << error["engy"] << " " << error["strs"] << " " << error["punish"]
       << endl;
    of.close();
    io.writeRadDist();
    io.writeAngDist();
    io.writeFrcDist();
    io.writeEngDist();
    for (int i = 0; i < 8; i++) {
      cout << "ebndl" << (dparams["fbndq"] = 0.03 + 0.003 * i) << endl;
      analyLoss();
    }
  }
  // data set
  // Melem aa;
  // for (int i = 0; i < 13; i++) {
  //   double delta = 0.02 * i;
  //   ItenT &tmp = aa.itdftm[delta]["opath"];
  //   cout << delta << " " << tmp.egy << " " << tmp.strm[1][1] << " "
  //        << tmp.strm[2][2] << " " << -0.1 * tmp.stsv[0] << " "
  //        << -0.1 * tmp.stsv[2] << " " << -0.1 * tmp.stsv[1] << " "
  //        << tmp.stsv[3] << " " << tmp.stsv[4] << " " << tmp.stsv[5] << endl;
  // }
  // check encoding
  // arma::mat v1 = encodev(mm);
  // arma::mat v2 = decodev(v1);
  // if (cmm.rank() == PFROOT)
  //   for (int i = 0; i < nvars; i++)
  //     cout << "i = " << i << " " << v1[i] << " " << v2[i] << " "
  //          << mm[i] - v2[i] << endl;
}

void pfHome::nloptGlobal(pfIO &io) {
  gcnt = 0;
  optdrv->prepOpt();
  optdrv->runOpt();
  io.writePot(ini);
}

void pfHome::doShift(pfForce &fcdrv, pfIO &io) {
  sparams["tmpfile"] = "dummy.tmp.0";
  io.writePot(ini);
  updaterhoMEAM(ini, fcdrv);
  double rho1 = oaverho;
  shiftRHO(ini);

  updaterhoMEAM(ini, fcdrv);
  double rho2 = oaverho;
  shiftEMF(rho2 - rho1);
  sparams["tmpfile"] = "dummy.tmp.1";
  io.writePot(ini);
}

double pfHome::errFunct(const vector<double> &x) {
  pfForce fcdrv(*this);
  double err = 0.0;
  int pfreq = 500;
  if (sparams["ptype"] == "EAM") {
    gcnt++;
    err = (fcdrv.*calobj[sparams["ptype"]])(x, 1);
    if (gcnt % pfreq == 0) printf("cnt %d err = %f \n", gcnt++, err);
  }
  return err;
}

double pfHome::errFunctGrad(const vector<double> &x, vector<double> &grad) {
  double y1, y2;
  vector<double> tmpx = x;
  for (unsigned int i = 0; i < x.size(); i++) {
    double dx = 1e-4 * x[i];
    tmpx[i] += dx;
    y1 = errFunct(tmpx);
    tmpx[i] -= (2 * dx);
    y2 = errFunct(tmpx);
    grad[i] = (y1 - y2) / (2 * dx);
  }
  return errFunct(x);
}
/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:31:59
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-03 16:22:23
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
#include "pfOptimizer.h"

using std::cout;
using std::endl;

void pfHome::increAnneal() {
  scnt = 0;
  int jobl[] = {PHI, PHI, RHO, MEAMF, MEAMG, PHI, RHO, PHI, RHO, MEAMF, MEAMG};
  FILE *fid;
  for (int i = 0; i < 10; i++) {
    simAnneal();
    upgrade(jobl[i]);
    iparams["kmax"]++;
    fid = fopen("record.txt", "a");
    fprintf(fid, "%d %f\n", i, recorderr[i]);
    fclose(fid);
  }
}

void pfHome::run(int argc, char *argv[]) {
  if (!sparams["opt"].compare("lmp")) {
    if (cmm.rank() == PFROOT) lmpdrv->calPhy();
  } else if (!sparams["opt"].compare("make") ||
             !sparams["opt"].compare("err")) {
    calErr();
    resample();
  } else if (!sparams["opt"].compare("gp"))
    GPsample();
  else if (!sparams["opt"].compare("anneal"))
    simAnneal();
  else if (!sparams["opt"].compare("evo"))
    diffEvo();
  else if (!sparams["opt"].compare("cmaes"))
    loopcmaes();
  else if (!sparams["opt"].compare("cnt"))
    cntcmaes();
  else if (!sparams["opt"].compare("incre"))
    increAnneal();
  else if (!sparams["opt"].compare("nlopt"))
    nloptGlobal();
  else if (!sparams["opt"].compare("shift"))
    doShift();
  else if (!sparams["opt"].compare("buildD03"))
    buildD03("d03", 7.400, 0.005);
  else if (!sparams["opt"].compare("anlz")) {
    calErr();
    calLat("bcc", 20);
    calLat("fcc", 20);
    lmpdrv->calLatticeBCC();
    exprs["lat"] = lmpdrv->exprs["lat"];
    calElas(25);

    lmpdrv->calLatticeFCC();
    lmpdrv->calLatticeHCP();
    lmpdrv->calSurfaceUrlx();
    lmpdrv->calGSFUrlx();

    remove("no");
    remove("log.lammps");
    remove("restart.equil");
    lmpdrv->exprs["bcc2hcp"] = lmpdrv->exprs["ehcp"] - lmpdrv->exprs["ebcc"];
    lmpdrv->exprs["bcc2fcc"] = lmpdrv->exprs["efcc"] - lmpdrv->exprs["ebcc"];

    vector<string> aa(
        {"lat", "bcc2fcc", "bcc2hcp", "suf110", "suf100", "suf111"});
    for (auto ee : aa) cout << ee << " " << lmpdrv->exprs[ee] << endl;

    for (int j : lmpdrv->gsfpnts)
      cout << std::setprecision(3) << lmpdrv->lgsf["111z110"][j] << " "
           << lmpdrv->lgsf["111e110"][j] << endl;
    for (int j : lmpdrv->gsfpnts)
      cout << std::setprecision(3) << lmpdrv->lgsf["111z211"][j] << " "
           << lmpdrv->lgsf["111e211"][j] << endl;
  }
}

void pfHome::calErr() {  // make potential
  arma::mat mm(nvars, 1, arma::fill::randu);
  for (int i = 0; i < nvars; i++) mm[i] = ini[i];
  (this->*calobj[sparams["ptype"]])(mm, 1);
  if (cmm.rank() == PFROOT) {
    double err = (this->*calobj[sparams["ptype"]])(mm, 1);
    (this->*calobj[sparams["ptype"]])(mm, EXT);
    (this->*write[sparams["ptype"]])();

    ofstream of("err.txt", std::ofstream::out);
    of << std::setprecision(4) << err << " " << error["frc"] << " "
       << error["engy"] << " " << error["strs"] << " " << error["punish"]
       << endl;
    of.close();
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
  // check the encoding
  // arma::mat v1 = encodev(mm);
  // arma::mat v2 = decodev(v1);
  // if (cmm.rank() == PFROOT)
  //   for (int i = 0; i < nvars; i++)
  //     cout << "i = " << i << " " << v1[i] << " " << v2[i] << " "
  //          << mm[i] - v2[i] << endl;
}

void pfHome::nloptGlobal() {
  gcnt = 0;
  optdrv->prepOpt();
  optdrv->runOpt();
  writePot(ini);
}

void pfHome::doShift() {
  cout << "before err = " << forceMEAMS(ini) << endl;

  sparams["tmpfile"] = "dummy.tmp.0";
  writePot(ini);

  updaterhoMEAM(ini);
  double rho1 = oaverho;
  shiftRHO(ini);

  updaterhoMEAM(ini);
  double rho2 = oaverho;

  shiftEMF(rho2 - rho1);

  sparams["tmpfile"] = "dummy.tmp.1";
  writePot(ini);

  cout << "after err = " << forceMEAMS(ini) << endl;
}

double pfHome::errFunct(const vector<double> &x) {
  double err = 0.0;
  int pfreq = 500;
  if (sparams["ptype"] == "EAM") {
    gcnt++;
    err = forceEAM(x);
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
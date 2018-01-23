/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:31:59
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-22 12:41:47
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
  if (!sparams["opt"].compare("lmp"))
    lmpdrv->calPhy();
  else if (!sparams["opt"].compare("make") || !sparams["opt"].compare("err"))
    calErr();  // buildD03("d03", 7.400, 0.005);
  else if (!sparams["opt"].compare("anneal"))
    simAnnealEAM();
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
  else if (!sparams["opt"].compare("anlz")) {
    calErr();
    readLmpMEAM();
    calLat("bcc");
    calElas();
    calPV();
    calSurf();
  }
}

void pfHome::calErr() {  // make potential
  arma::mat mm(nvars, 1, arma::fill::randu);
  for (int i = 0; i < nvars; i++) mm[i] = ini[i];
  if (!sparams["ptype"].compare("EAM")) {
    forceEAM(mm, 1);
    if (cmm.rank() == PFROOT) {
      for (int i = 0; i < 10; i++) {
        double err = forceEAM(mm, 1);
        cout << "Err " << i << " " << err << endl;
      }
      forceEAM(mm, EXT);
    }
  } else if (!sparams["ptype"].compare("ADP")) {
    cout << "error = " << forceADP(ini, 0) << endl;
  } else if (!sparams["ptype"].compare("MEAM"))
    cout << "error = " << forceMEAM(mm) << endl;
  // if (cmm.rank() == PFROOT)
  //   sparams["ptype"] == "MEAM" ? writeMEAM() : writeLMPS();
}

void pfHome::nloptGlobal() {
  gcnt = 0;
  optdrv->prepOpt();
  optdrv->runOpt();
  writePot(ini);
}

void pfHome::doShift() {
  cout << "before err = " << forceMEAM(ini) << endl;

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

  cout << "after err = " << forceMEAM(ini) << endl;
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
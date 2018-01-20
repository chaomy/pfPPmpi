/*
 * @Author: chaomy
 * @Date:   2018-01-15 00:24:43
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-16 14:00:50
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
#include "pfOptimizer.h"

pfHome::pfHome(int argc, char* argv[]) {
  ricut = 2.08;
  rocut = 5.05;
  ftn = tln = 0;
  for (int ii : {0, 1, 2}) mfrc[ii] = 0.0;
  iparams["atomicNum"] = 40;
  dparams["mass"] = 92.906400;
  sparams["tmpdir"] = string("dirtmp");
  sparams["lmpdir"] = string("dirlmp");
  outMkdir(sparams["tmpdir"]);
  outMkdir(sparams["lmpdir"]);

  sparams["tmpfile"] = string("dummy.tmp");
  sparams["parfile"] = string("dummy.param");
  sparams["cnffile"] = string("dummy.config");
  sparams["potfile"] = string("dummy.pot");
  sparams["lmppot"] = string("lmp.pot");

  // set the boundary
  double hilim[6] = {2., 0.5, -5., 0.5, 0.5, 20};
  double lolim[6] = {-0.5, -0.1, -10., -0.5, -0.5, -20};
  for (int i = 0; i < 6; i++) {
    gradRight.push_back(0);
    hil.push_back(hilim[i]);
    lol.push_back(lolim[i]);
  }

  lorho = 0.4;
  hirho = 1.0;

  gradRight[EMF] = 1;
  gradRight[MEAMG] = 1;
  parseArgs(argc, argv);
  pfInit();
  lmpdrv = new pfLMPdrv(argc, argv, this);
  optdrv = new pfOptimizer(this);
};

pfHome::~pfHome() {
  funcs.clear();
  ini.clear();
  lob.clear();
  hib.clear();
  configs.clear();
  startps.clear();
  if (lmpdrv) delete lmpdrv;
  if (optdrv) delete optdrv;
}

void pfHome::pfInit() {
  initParam();
  initTargs();
  readConfig();
  readPot();
  if (!sparams["ptype"].compare("MEAM")) {
    initNeighsFull();
    initAngles();
  } else
    initNeighs();
  forceDis();
}
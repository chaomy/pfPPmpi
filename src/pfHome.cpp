/*
 * @Author: chaomy
 * @Date:   2018-01-15 00:24:43
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-04-04 13:31:44
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
#include "pfOptimizer.h"
namespace mpi = boost::mpi;
#define MXEL 5

pfHome::pfHome(int argc, char* argv[])
    : ricut(2.08),
      rocut(6.00),
      mfrc(3),
      hil({5.0, 2.0, 0.0, 1.0, 1.0}),
      lol({-1.0, -30.0, -10.0, -0.5, -0.5}),
      nelt(1),
      Ec_meam(MXEL, vector<double>(MXEL, 1.55)),
      re_meam(MXEL, vector<double>(MXEL, 3.2)),
      Omega_meam(MXEL),
      Z_meam(MXEL),
      A_meam(MXEL),
      alpha_meam(MXEL, vector<double>(MXEL)),
      rho0_meam(MXEL),
      delta_meam(MXEL, vector<double>(MXEL, 0.0)),
      beta0_meam(MXEL),
      beta1_meam(MXEL),
      beta2_meam(MXEL),
      beta3_meam(MXEL),
      t0_meam(MXEL),
      t1_meam(MXEL),
      t2_meam(MXEL),
      t3_meam(MXEL),
      rho_ref_meam(MXEL),
      ibar_meam(MXEL),
      ielt_meam(MXEL),
      lattce_meam(MXEL, vector<lattice_t>(MXEL)),
      nn2_meam(MXEL, vector<int>(MXEL, 1)),
      zbl_meam(MXEL, vector<int>(MXEL, 0)),
      eltind(MXEL, vector<int>(MXEL)),
      attrac_meam(MXEL, vector<double>(MXEL, 0.0)),
      repuls_meam(MXEL, vector<double>(MXEL, 0.0)),
      Cmin_meam(MXEL, vector<vector<double>>(MXEL, vector<double>(MXEL, 0.49))),
      Cmax_meam(MXEL, vector<vector<double>>(MXEL, vector<double>(MXEL, 2.8))),
      ebound_meam(MXEL,
                  vector<double>(MXEL, pow(2.8, 2) / (4.0 * (2.8 - 1.0)))),
      rc_meam(4.8),
      delr_meam(0.1),
      gsmooth_factor(99.0),
      augt1(0),
      ialloy(1),
      mix_ref_t(0),
      emb_lin_neg(0),
      bkgd_dyn(0),
      erose_form(2) {
  calfrc["EAM"] = &pfHome::forceEAM;
  calfrc["MEAMS"] = &pfHome::forceMEAMS;
  calfrc["MEAMC"] = &pfHome::forceMEAMC;

  calobj["EAM"] = &pfHome::forceEAM;
  calobj["MEAMS"] = &pfHome::forceMEAMS;
  calobj["MEAMC"] = &pfHome::forceMEAMC;

  write["EAM"] = &pfHome::writeLMPS;
  write["TMP"] = &pfHome::writePot;
  write["MEAMS"] = &pfHome::writeMEAMS;
  write["MEAMC"] = &pfHome::writeMEAMC;

  read["EAMS"] = &pfHome::readPot;
  read["MEAMS"] = &pfHome::readMEAMS;
  read["MEAMC"] = &pfHome::readMEAMC;

  build["bcc"] = &pfHome::buildbccPrim;
  build["fcc"] = &pfHome::buildfccPrim;

  latticemp = vector<string>(
      {"fcc", "bcc", "hcp", "dim", "dia", "b1", "c11", "l12", "b2"});

  if (cmm.rank() == PFROOT) {
    ftn = tln = 0;
    for (int ii : {0, 1, 2}) mfrc[ii] = 0.0;
    iparams["atomicNum"] = 40;
    dparams["mass"] = 92.906400;
    outMkdir(sparams["lmpdir"] = string("dirlmp"));
    // sparams["tmpdir"] = string("dirtmp");
    // outMkdir(sparams["tmpdir"]);

    sparams["tmpfile"] = string("pf.tmp");
    sparams["parfile"] = string("pf.par");
    sparams["cnffile"] = string("pf.cnf");
    sparams["lmppot"] = string("pf.lmp");
    sparams["potfile"] = string("meam.lib");
    sparams["meamcnt"] = string("meam.cnt");
    sparams["meamlib"] = string("meam.tmp");
    sparams["meampar"] = string("meam.param");

    lorho = 0.4, hirho = 1.0;
    parseArgs(argc, argv);
    pfInit();
  }

  (cmm.barrier)();  //  important, be careful
  bcdata();
  if (!sparams["ptype"].compare("MEAMS")) {
    initNeighsFull();
    initAngles();
  } else if (!sparams["ptype"].compare("MEAMC")) {
    initNeighsFull();
  } else
    initNeighs();

  cmmlm = cmm.split(cmm.rank() == PFROOT);  // split group lm to run lammps
  lmpdrv = new pfLMPdrv(argc, argv, this);
  assignConfigs(2);
  (cmm.barrier)();  //  important!

  if (cmm.rank() == 2)
    for (auto cc : configs) {
      cout << "cnf weigh " << cc.weigh << endl;
      if (cc.weigh > 100) exit(1);
    }
  // temporarily close these functionalities
  // optdrv = new pfOptimizer(this);
};

pfHome::~pfHome() {
  funcs.clear();
  ini.clear();
  lob.clear();
  hib.clear();
  configs.clear();
  startps.clear();
  if (lmpdrv) delete lmpdrv;
  // if (optdrv) delete optdrv;
}

void pfHome::assignConfigs(int tag) {
  if (tag == 1) {  // evenly assign configurations
    int del = nconfs / cmm.size();
    locstt = del * cmm.rank();
    locend = del * (cmm.rank() + 1);
    for (int i = 0; i < nconfs; i++) {
      if (i % cmm.size() == cmm.rank()) locls.push_back(i);
    }
  }
  if (tag == 2) {  // assign configurations based on number of atoms
    locnatoms = vector<int>(cmm.size(), 0);
    for (int i = 0; i < nconfs; i++) {
      int mid = 0;
      int mnatms = INT_MAX;
      for (int j = 0; j < cmm.size(); j++) {
        if (locnatoms[j] <= mnatms) {
          mnatms = locnatoms[j];
          mid = j;
        }
      }
      locnatoms[mid] += configs[i].atoms.size();
      if (cmm.rank() == mid) locls.push_back(i);
    }
  }
  cout << " rank " << cmm.rank() << " has " << locnatoms[cmm.rank()] << " atoms"
       << endl;
}

void pfHome::pfInit() {
  initParam();
  initTargs();
  readConfig();
  (this->*read[sparams["ptype"]])();
}
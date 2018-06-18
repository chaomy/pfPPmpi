/*
 * @Author: chaomy
 * @Date:   2018-01-15 00:24:43
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-18 16:40:28
 */

#include "pfConf.h"
#include "pfForce.h"
#include "pfIO.h"
#include "pfLmpDrv.h"
#include "pfMEAMC.h"
#include "pfOptimizer.h"
#include "pfPhy.h"

namespace mpi = boost::mpi;

pfHome::pfHome(int argc, char* argv[])
    : ricut(2.08),
      rocut(6.00),
      mfrc(3),
      hil({5.0, 2.0, 0.0, 1.0, 1.0}),
      lol({-1.0, -30.0, -10.0, -0.5, -0.5}) {
  calfrc["EAM"] = &pfHome::pfForce::forceEAM;
  calfrc["MEAMS"] = &pfHome::pfForce::forceMEAMS;
  // calfrc["MEAMC"] = &pfHome::pfForce::pfMEAMC::forceMEAMC;

  calobj["EAM"] = &pfHome::pfForce::forceEAM;
  calobj["MEAMS"] = &pfHome::pfForce::forceMEAMS;
  // calobj["MEAMC"] = &pfHome::pfForce::pfMEAMC::forceMEAMC;

  write["EAM"] = &pfHome::pfIO::writeLMPS;
  write["TMP"] = &pfHome::pfIO::writePot;
  write["MEAMS"] = &pfHome::pfIO::writeMEAMS;

  read["EAMS"] = &pfHome::pfIO::readPot;
  read["MEAMS"] = &pfHome::pfIO::readMEAMS;

  cout << "hello ?" << endl;

  pfPhy phdrv(*this);
  pfConf cdrv(*this);
  pfForce fcdrv(*this);
  pfIO io(*this);

  if (cmm.rank() == PFROOT) {
    ftn = tln = 0;
    for (int ii : {0, 1, 2}) mfrc[ii] = 0.0;
    iparams["atomicNum"] = 40;
    dparams["mass"] = 92.906400;
    io.outMkdir(sparams["lmpdir"] = string("dirlmp"));
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
    pfInit(io);
  }

  (cmm.barrier)();
  bcdata();
  if (!sparams["ptype"].compare("MEAMS")) {
    cdrv.initNeighsFull();
    cdrv.initAngles();
  } else if (!sparams["ptype"].compare("MEAMC")) {
    cdrv.initNeighsFull();
  } else
    cdrv.initNeighs();

  cmmlm = cmm.split(cmm.rank() == PFROOT);  // split group lm to run lammps
  lmpdrv = new pfLMPdrv(argc, argv, *this);
  assignConfigs(2);
  (cmm.barrier)();  //  important!
  // optdrv = new pfOptimizer(this); // temporarily deactivate functionalities

  run(argc, argv, fcdrv, cdrv, phdrv, io);
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
    for (int i = 0; i < nconfs; i++)
      if (i % cmm.size() == cmm.rank()) locls.push_back(i);
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

void pfHome::pfInit(pfIO& io) {
  initParam(io);
  initTargs();
  io.readConfig();
  (io.*read[sparams["ptype"]])();
}
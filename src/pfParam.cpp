/*
 * @Author: chaomy
 * @Date:   2017-11-05 22:29:46
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-16 16:36:29
 */

#include "pfIO.h"

void pfHome::initParam(pfIO& io) {
  sparams["spline"] = string("nat");
  sparams["elem"] = string("Nb");
  sparams["ptype"] = string("MEAMC");
  sparams["pairstyle"] = string("meam/spline");
  sparams["alg"] = string("LN_SBPLX");
  sparams["opt"] = string("nlopt");

  dparams["temp"] = 0.01;
  dparams["istep"] = 0.3;
  dparams["ivari"] = 0.3;
  dparams["pratio"] = 0.5;
  dparams["eweight"] = 5.0;
  dparams["sweight"] = 1.0;
  dparams["pweight"] = 1.0;
  dparams["pshift"] = 10;
  dparams["xtol"] = 1e-4;
  dparams["ftol"] = 1e-6;
  dparams["bwidth"] = 10.0;
  dparams["fbndl"] = 30 * (dparams["fbndq"] = 0.02);
  dparams["ebndl"] = 30 * (dparams["ebndq"] = 0.01);
  dparams["frceps"] = 0.1;

  iparams["maxstep"] = 10000;
  iparams["resfreq"] = 10;
  iparams["lmpfreq"] = 50;
  iparams["kmax"] = 1000;  // number of outer loop in simulated annealing
  iparams["runlmp"] = 0;
  io.readParam();
  sparams["lmpfile"] = string("lmp.") + sparams["ptype"];
}

void pfHome::pfIO::readParam() {
  ifstream fid;
  fid.open(sparams["parfile"].c_str());
  if (!fid.is_open()) cerr << " error opening " + sparams["parfile"] << endl;
  vector<string> segs(1, " ");
  string buff;

  while (getline(fid, buff)) {
    segs.clear();
    cout << segs[0] << " " << segs[1] << endl;
    if (!segs[0].compare("alg"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("elem"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("opt"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("ptype"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("pairstyle"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("spline"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("maxstep"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("resfreq"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("lmpfreq"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("runlmp"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("kmax"))
      iparams[segs[0]] = stoi(segs[1]);
    else
      dparams[segs[0]] = stof(segs[1]);
  }
  fid.close();
}

/**** parse the command line ****/
void pfHome::parseArgs(int argc, char* argv[]) {
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      sparams["parfile"] = string(argv[++i]);
    else if (!strcmp(argv[i], "--c") || !strcmp(argv[i], "-c"))
      sparams["cnffile"] = string(argv[++i]);
    else if (!strcmp(argv[i], "--f") || !strcmp(argv[i], "-f"))
      sparams["potfile"] = string(argv[++i]);
  }
}
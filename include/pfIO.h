#ifndef PF_IO_H
#define PF_IO_H

#include "pfHome.h"

class pfHome::pfIO {
  pfHome& hm;
  unordered_map<string, double>& dparams;  // double parameters
  unordered_map<string, int>& iparams;     // integer parameters
  unordered_map<string, string>& sparams;  // string parameters
  vector<Config>& configs;                 // configurations
  vector<Func>& funcs;                     // spline functions
  vector<double>& ini;                     // hold spline node values
  vector<int>& optidx;                     // functions that to be optimized
  vector<int>& smthidx;
  vector<double>& hil;   // 5 + 1  default max values of each function
  vector<double>& lol;   // 5 + 1  default min values of each function
  vector<int>& startps;  // used in annealing method to update functions
  vector<int>& locls;
  int& tln;  // total number of atoms
  int& nconfs;
  int& nfuncs;

 public:
  pfIO(pfHome& x)
      : hm(x),
        dparams(x.dparams),
        iparams(x.iparams),
        sparams(x.sparams),
        configs(x.configs),
        funcs(x.funcs),
        ini(x.ini),
        optidx(x.optidx),
        smthidx(x.smthidx),
        hil(x.hil),
        lol(x.lol),
        startps(x.startps),
        locls(x.locls),
        tln(x.tln),
        nconfs(x.nconfs),
        nfuncs(x.nfuncs) {}

  // inputs
  void readConfig();
  void readPOSCAR();
  void readPot();
  void readMEAMC();
  void readMEAMS();
  void readMEAMCcnt();
  void readParam();

  // outputs
  void writePot();
  void writePot(const vector<double>& vv);
  void writePot(const string& s);
  void writeLMPS();
  void writeLMPS(const vector<double>& vv);
  void writePOSCAR(const Config& cc, string fnm = "POSCAR.vasp");
  void writeMEAMS();
  void writeRadDist();
  void writeAngDist();
  void writeFrcDist();
  void writeEngDist();

  // utils
  void outMkdir(string mdir);
};

#endif  // PF_IO_H

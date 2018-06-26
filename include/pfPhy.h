#ifndef PF_PHY_H
#define PF_PHY_H

#include "pfHome.h"

class pfHome::pfPhy {
  pfHome& hm;
  Melem& mele;
  unordered_map<string, string>& sparams;  // string parameters
  unordered_map<string, double>& pwgh;     // pwght of errors
  unordered_map<string, double>& exprs;    // lammps errors
  unordered_map<string, double>& targs;    // target values
  unordered_map<string, double>& perr;     // physical errors
  Ctrib& error;
  unordered_map<string, void (pfHome::pfForce::*)(Config&)>& calfrc;

 public:
  pfPhy(pfHome& x)
      : hm(x),
        mele(x.mele),
        sparams(x.sparams),
        pwgh(x.pwgh),
        exprs(x.exprs),
        targs(x.targs),
        perr(x.perr),
        error(x.error),
        calfrc(x.calfrc) {}

  void calLat(string key, pfForce& f, pfConf& c);
  void calLat(string key, int n, pfForce& f, pfConf& c);
  void calPV(pfForce& f, pfConf& c);
  void calElas(pfForce& f, pfConf& c);
  void calElas(int npts, pfForce& f, pfConf& c);
  void calSurf(pfForce& f, pfConf& c);
};

#endif  // PF_PHY_H
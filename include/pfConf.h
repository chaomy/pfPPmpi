#ifndef PF_CONF_H
#define PF_CONF_H

#include "pfHome.h"

class pfHome::pfConf {
 private:
  /* to examine potential */
  Config ubcc;  // primitive bcc
  Config cbcc;  // conventional bcc
  unordered_map<string, vector<Config>> mpcf;
  unordered_map<string, vector<double>> mpvc;
  unordered_map<string, Config (pfHome::pfConf::*)(const double& lat)> build;

  /**/
  pfHome& hm;
  mpi::communicator& cmm;
  vector<Config>& configs;                 // configurations
  vector<Func>& funcs;                     // spline functions
  unordered_map<string, string>& sparams;  // string parameters
  double& ricut;
  double& rocut;
  int& nconfs;
  int& ftn;

 public:
  pfConf(pfHome& x)
      : hm(x),
        cmm(x.cmm),
        configs(x.configs),
        funcs(x.funcs),
        sparams(x.sparams),
        ricut(x.ricut),
        rocut(x.rocut),
        nconfs(x.nconfs),
        ftn(x.ftn) {
    build["bcc"] = &pfHome::pfConf::buildbccPrim;
    build["fcc"] = &pfHome::pfConf::buildfccPrim;
  }

 public:
  // initialize box
  void initBox(Config& cc);

  void wrapAtomPos(Config& cc);

  // initialize neighbors
  void initNeighs();
  void initNeighs(Config& cc);
  void initNeighsFull();           /* meam */
  void initNeighsFull(Config& cc); /* meam */

  // initialize angles
  void initAngles();
  void initAngles(Config& cc);
  void initAnglesSameCutOff();

  // set neighbors
  void setNeighslot(Neigh& n, Func f, double r);
  void setNeighslotStd(Neigh& n, Func f, double r);
  void updateNeighslot(Neigh& n, Func f, double r, int id);

  // set angles
  void setAngleslot(Angle& a, Func f, double r);
  void setAngleslotStd(Angle& a, Func f, double r);

  Config addvolm(const Config&, const double& dl);
  Config addotho(const Config&, const double& dl);
  Config addmono(const Config&, const double& dl);
  Config addstrain(Config cc, const vector<vector<double>>& str);

  // build configs
  void buildbcc(const string& kk, const double& gs, const double& dl);
  void buildfcc(const string& kk, const double& gs, const double& dl);
  void buildhcp(const string& kk, const double& gs, const double& dl);
  void buildD03(const string& kk, const double& gs, const double& dl);

  Config buildbccConv(const double& lat);
  Config buildfccConv(const double& lat);
  Config buildbccPrim(const double& lat);
  Config buildfccPrim(const double& lat);
  Config buildsur100(const double& lat, const string& tag);
  Config buildsur110(const double& lat, const string& tag);
  Config buildsur211(const double& lat, const string& tag);
  Config buildhcp(const double& la, const double& lc);
  Config buildD03(const double& lat);

  friend class pfHome::pfPhy;
};

#endif  // PF_CONF_H

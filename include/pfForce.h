#ifndef pfForce_H_
#define pfForce_H_

#include "pfHome.h"

class pfHome::pfForce {
 private:
  pfHome& hm;
  vector<Func>& funcs;
  vector<Config>& configs;
  mpi::communicator& cmm;
  unordered_map<string, double>& dparams;
  vector<double>& ini;
  vector<int>& smthidx;
  vector<int>& optidx;
  vector<int>& locls;
  Ctrib& error;
  Ctrib& weigh;
  int& locstt;
  int& locend;
  double& omaxrho;
  double& ominrho;

 public:
  class pfMEAMC;
  pfForce(pfHome& x)
      : hm(x),
        funcs(x.funcs),
        configs(x.configs),
        cmm(x.cmm),
        dparams(x.dparams),
        ini(x.ini),
        smthidx(x.smthidx),
        optidx(x.optidx),
        locls(x.locls),
        error(x.error),
        weigh(x.weigh),
        locstt(x.locstt),
        locend(x.locend),
        omaxrho(x.omaxrho),
        ominrho(x.ominrho) {}
  double forceEAM(vector<Func>& ffs, int tag);
  double forceEAM(const arma::mat& vv);
  double forceEAM(const arma::mat& vv, int tg);
  double forceADP(const vector<double>& vv, int tag);
  double forceMEAMS(const arma::mat& vv);
  double forceMEAMS(const arma::mat& vv, int tg);
  double forceMEAMSNoForce(const arma::mat& vv, int tg);
  double forceMEAMSStress(const arma::mat& vv, int tg);
  double forceMEAMSStressPunish(const arma::mat& vv, int tg);
  // double forceMEAMC(const arma::mat& vv, int tg);
  // void forceMEAMC(Config& cc);
  void forceMEAMS(Config& cc);
  void forceMEAMSNoForce(Config& cc);
  void forceMEAMSStress(Config& cc);
  void forceEAM(Config& cc);
  void stressMEAM(Config& cc);

  // spline interpolation
  void spltra(Func& func, double r, double& val, double& grad);
  void spltrai(Func& func, double r, double& val, double& grad);
  void splineNe(Func& func, int flag);
  void splineEd(Func& func, int flag);

  void splintEd(const Func& func, double r, double& val);  // no gradients
  void splintEd(const Func& func, double r, double& val, double& grad);

  void splint(const Func& func, double r, double& val);  // no gradients
  void splint(const Func& func, double r, double& val, double& grad);

  void splint(const Func& func, int k, double b, double step,
              double& val);  // no gradients
  void splint(const Func& func, int k, double b, double step, double& val,
              double& grad);
};
#endif  // pfForce_H_
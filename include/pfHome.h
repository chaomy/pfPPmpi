#ifndef pfHome_H_
#define pfHome_H_

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include "armadillo"
#include "lmpMatrix.h"
#include "math_special.h"
#include "nlopt.hpp"
#include "pfDefines.h"
#include "pfEle.h"
#include "spline.h"

using std::abs;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::to_string;
using std::unordered_map;
using std::vector;

namespace mpi = boost::mpi;

class pfOptimizer;
class pfHome {
 public:
  mpi::environment env;
  mpi::communicator cmm;
  mpi::communicator cmmlm;

  class pfLMPdrv;
  class pfForce;
  class pfPhy;
  class pfConf;
  class pfIO;

 private:
  int ftn;   // number of atoms used for fitting
  int tln;   // total number of atoms
  int gcnt;  // count the times of calling force calculation
  int scnt;  // count the times of calling simulated annealing
  int nvars;
  int nfuncs;
  int nconfs;
  int locstt, locend;
  vector<int> locnatoms;
  vector<int> locls;
  vector<int> smthidx;
  double fsm, phyweigh;

  double ricut;
  double rocut;
  double rhcut;
  double lorho;
  double hirho;
  double punish;
  double physic;

  double ominrho;  // overall min density
  double omaxrho;  // overall max density
  double oaverho;  // overall average density

  /* parameters */
  unordered_map<string, double> dparams;            // double parameters
  unordered_map<string, int> iparams;               // integer parameters
  unordered_map<string, string> sparams;            // string parameters
  unordered_map<string, vector<double>> meamparms;  // hold MEAMC parameters
  vector<string> elems;                             //  element name

  unordered_map<string, double> targs;  // target values
  unordered_map<string, double> exprs;  // lammps errors
  unordered_map<string, double> weigh;  // weight of errors
  unordered_map<string, double> error;  // errors

  /* mapping functions */
  unordered_map<string, void (pfHome::pfForce::*)(Config&)> calfrc;
  unordered_map<string,
                double (pfHome::pfForce::*)(const arma::mat& vv, int tg)>
      calobj;
  unordered_map<string, void (pfHome::pfIO::*)()> write;
  unordered_map<string, void (pfHome::pfIO::*)()> read;

  vector<Config> configs;  // configurations
  vector<Func> funcs;      // spline functions
  // vector<Func> fprec;  to be deleted

  vector<int> startps;  // used in annealing method to update functions
  vector<int> endps;    // used in annealing method to update functions
  pfLMPdrv* lmpdrv;     // lammps driver
  pfOptimizer* optdrv;  // optimizer
  Melem mele;           // element data

  vector<double> mfrc;       // min force of each force vectors
  vector<double> hil;        // 5 + 1  default max values of each function
  vector<double> lol;        // 5 + 1  default min values of each function
  vector<double> recorderr;  // growing variable, record error of each run
  vector<double> ini;        // hold spline node values
  vector<double> hib;        // hold high bound spline node values
  vector<double> lob;        // hold low bound spline node values
  vector<double> deb;        // hib - lob
  vector<int> optidx;        // functions that to be optimized

 public:
  pfHome(int argc, char* argv[]);
  ~pfHome();

  unordered_map<string, string> gsparams() const { return sparams; };
  unordered_map<string, double> gdparams() const { return dparams; };
  unordered_map<string, int> giparams() const { return iparams; };

  // initialization
  void pfInit(pfIO& io);
  void parseArgs(int argc, char* argv[]);
  void initParam(pfIO& io);
  void initTargs();
  void assignConfigs(int tag = 2);
  void setSplineVariables();

  // force calculation
  void doShift(pfForce& f, pfIO& i);  // make the rho / emf good
  void run(int argc, char* argv[], pfForce&, pfConf&, pfPhy&, pfIO&);
  void calErr(pfForce& f, pfIO& i);
  void calErr(int tm);  // for debugging
  void updaterho(vector<double>& vv);
  void updaterhoMEAM(vector<double>& vv, pfForce& fcdrv);
  void resample(pfIO&);
  double errFunct(const vector<double>& x);
  double errFunctGrad(const vector<double>& x, vector<double>& g);

  // optimization
  double encode(const double& val, const int& idx);
  double decode(const double& val, const int& idx);
  arma::mat encodev(const arma::mat& vv);
  arma::mat encodev(const vector<double>& vv);
  arma::mat decodev(const arma::mat& vv);
  vector<double> decodestdv(const arma::mat& vv);
  vector<double> decodestdv(const vector<double>& vv);
  vector<double> encodestdv(const vector<double>& vv);
  void simAnneal(pfForce& f, pfIO& io);
  void simAnnealSpline(pfForce& f, pfIO& io);
  void randomize(vector<double>& vv, const int n, const vector<double>& v);
  void randomizeSpline(vector<double>& vv, const int n,
                       const vector<double>& v);
  int rescaleEMF(vector<double>& vv, pfForce& fcdrv);
  int rescaleEMF(arma::mat& vv);
  int rescaleRHO(vector<double>& vv);
  bool checkBoundary(const arma::mat& vv);
  void checkupdateBoundary(arma::mat& vv);
  void updateBoundary(const arma::mat& vv);
  void shiftRHO(vector<double>& vv);
  void shiftEMF(double shift);
  void nloptGlobal(pfIO&);

  // Diff evo
  void initPop();
  void diffEvo(pfForce& f, pfIO& io);

  // Gaussian Process
  void GPsample(pfIO&);

  // CMA-ES
  void loopcmaes(pfForce& f, pfIO& io);
  void cntcmaes(pfForce& f, pfIO& io);
  double cmaes(arma::mat& iterate, pfForce& f, pfIO& io);
  double testFunc(arma::mat& coordinates);

  // random
  double randNormal();
  double randUniform();
  double randUniform(const double min, const double max);

  void upgrade(int id, pfForce& f, pfConf& c);       // increase nodes
  void increAnneal(pfForce& f, pfConf& c, pfIO& i);  // increasing nodes
  void recordStage(int cnt, pfIO& io);

  // MPI utilitis
  void bcdata(); /* broadcast data and parameters */

  // use it to adjust parameters in fitting
  void deleteAtoms();
  void cutoffNeighs();
  void loopBwth();
  void forceDis();
  void analyLoss();
  void lmpCheck(int i, ofstream& of1);

  friend class pfOptimizer;
};

void inline split(const string& s, const char* delim, vector<string>& v) {
  // duplicate original string, return a char pointer and free  memories
  char* dup = strdup(s.c_str());
  char* token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline arma::mat pfHome::decodev(const arma::mat& vv) {
  arma::mat rs(nvars, 1);
  for (int i = 0; i < nvars; i++)
    rs[i] = lob[i] + deb[i] * 0.5 * (1. - cos(PI * 0.1 * vv[i]));
  return rs;
}

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline double pfHome::decode(const double& val, const int& i) {
  return lob[i] + deb[i] * 0.5 * (1. - cos(PI * 0.1 * val));
}

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline vector<double> pfHome::decodestdv(const vector<double>& vv) {
  vector<double> rs(vv.size());
  for (int i = 0; i < nvars; i++)
    rs[i] = lob[i] + deb[i] * 0.5 * (1. - cos(PI * 0.1 * vv[i]));
  return rs;
}

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline vector<double> pfHome::decodestdv(const arma::mat& vv) {
  vector<double> rs(vv.size());
  for (int i = 0; i < nvars; i++)
    rs[i] = lob[i] + deb[i] * 0.5 * (1. - cos(PI * 0.1 * vv[i]));
  return rs;
}

// [a, b] -> [0, 10]
inline arma::mat pfHome::encodev(const arma::mat& vv) {
  arma::mat rs(nvars, 1);
  for (int i = 0; i < nvars; i++)
    rs[i] = 10 * acos(1. - 2. / deb[i] * (vv[i] - lob[i])) * INVPI;
  return rs;
}

inline arma::mat pfHome::encodev(const vector<double>& vv) {
  arma::mat rs(nvars, 1);
  for (int i = 0; i < nvars; i++)
    rs[i] = 10 * acos(1. - 2. / deb[i] * (vv[i] - lob[i])) * INVPI;
  return rs;
}

inline double pfHome::encode(const double& val, const int& i) {
  return 10 * acos(1. - 2. / deb[i] * (val - lob[i])) * INVPI;
}

// [a, b] -> [0, 10]
inline vector<double> pfHome::encodestdv(const vector<double>& vv) {
  vector<double> rs(vv.size());
  for (int i = 0; i < nvars; i++)
    rs[i] = 10 * acos(1. - 2. / deb[i] * (vv[i] - lob[i])) * INVPI;
  return rs;
}

template <typename TYPE>
static inline void setall2d(vector<vector<TYPE>>& vv, const TYPE val) {
  for (auto& v1 : vv)
    for (auto& v2 : v1) v2 = val;
}

template <typename TYPE>
static inline void setall3d(vector<vector<vector<TYPE>>> vv, const TYPE val) {
  for (auto& v1 : vv)
    for (auto& v2 : v1)
      for (auto& v3 : v2) v3 = val;
}

inline bool pfHome::checkBoundary(const arma::mat& vv) {
  for (int i = 0; i < nvars; i++)
    if (vv[i] < 0.1 || vv[i] > 9.9) return true;
  return false;
}

inline void pfHome::checkupdateBoundary(arma::mat& vv) {
  for (int i = 0; i < nvars; i++) {
    if (vv[i] < 0.1 || vv[i] > 9.9) {
      double val = decode(vv[i], i);
      double vari = (fabs(val) >= 1e-8) ? dparams["ivari"] * fabs(val) : 0.001;
      cout << i << " before update boundary" << lob[i] << " " << hib[i] << endl;
      lob[i] = val - vari;
      hib[i] = val + vari;
      deb[i] = 2 * vari;
      cout << i << " after update boundary" << lob[i] << " " << hib[i] << endl;
      vv[i] = 5.0;
    }
  }
}

inline void pfHome::updateBoundary(const arma::mat& vv) {
  for (int i = 0; i < nvars; i++) {
    double vari =
        (fabs(vv[i]) >= 1e-8) ? dparams["ivari"] * fabs(vv[i]) : 0.001;
    lob[i] = vv[i] - vari;
    hib[i] = vv[i] + vari;
    deb[i] = 2 * vari;
  }
}

#endif  // pfHome_H_
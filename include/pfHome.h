#ifndef pfHome_H_
#define pfHome_H_

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include "armadillo"
#include "lmpMatrix.h"
#include "nlopt.hpp"
#include "pfDefines.h"
#include "pfEle.h"
#include "spline.h"

using std::string;
using std::unordered_map;
using std::vector;
using namespace MMatrix;

class pfLMPdrv;
class pfOptimizer;

class pfHome {
 private:
  int ftn;  // number of atoms used for fitting
  int tln;  // total number of atoms
  int chid;
  int gcnt;  // count the times of calling force calculation
  int scnt;  // count the times of calling simulated annealing
  int nvars;
  int nfuncs;
  int nconfs;
  int locstt;
  int locend;
  double mfrc[3];
  double fsm;

  double ricut;
  double rocut;
  double rhcut;
  double lorho;
  double hirho;
  double punish;
  double physic;

  double ominrho;
  double omaxrho;
  double oaverho;

  /* parameters */
  unordered_map<string, double> dparams;
  unordered_map<string, int> iparams;
  unordered_map<string, string> sparams;

  /* tests */
  Config ubcc;  // primitive bcc
  Config cbcc;  // conventional bcc
  unordered_map<string, vector<Config>> mpcf;
  unordered_map<string, vector<double>> mpvc;

  unordered_map<string, double> targs;
  unordered_map<string, double> exprs;
  unordered_map<string, double> weigh;
  unordered_map<string, double> error;

  vector<Config> configs;
  vector<Func> funcs;
  vector<Func> fprec;

  vector<int> startps;
  vector<int> endps;
  vector<int> gradRight;
  pfLMPdrv* lmpdrv;
  pfOptimizer* optdrv;
  Melem mele;

  vector<double> hil;  // 5 + 1
  vector<double> lol;  // 5 + 1
  vector<double> recorderr;
  vector<double> ini;
  vector<double> lob;
  vector<double> hib;

 public:
  pfHome(int argc, char* argv[]);
  ~pfHome();

  unordered_map<string, string> gsparams() const { return sparams; };
  unordered_map<string, double> gdparams() const { return dparams; };
  unordered_map<string, int> giparams() const { return iparams; };

  // initialization
  void pfInit();
  void parseArgs(int argc, char* argv[]);
  void initParam();
  void initTargs();

  void initBox(Config& cc);
  void initNeighs();
  void initNeighs(Config& cc);
  void initNeighsFull();           /* meam */
  void initNeighsFull(Config& cc); /* meam */
  void initAngles();
  void initAngles(Config& cc);
  void initAnglesSameCutOff();
  void setNeighslot(Neigh& n, Func f, double r);
  void setNeighslotStd(Neigh& n, Func f, double r);
  void updateNeighslot(Neigh& n, Func f, double r, int id);
  void setAngleslot(Angle& a, Func f, double r);
  void setAngleslotStd(Angle& a, Func f, double r);

  // spline interpolation
  void spltra(Func& func, double r, double& val, double& grad);
  void spltrai(Func& func, double r, double& val, double& grad);
  void splineNe(Func& func, int flag);
  void splineEd(Func& func, int flag);
  void splintEd(const Func& func, double r, double& val);
  void splintEd(const Func& func, double r, double& val, double& grad);
  void splint(const Func& func, double r, double& val, double& grad);
  void splint(const Func& func, double r, double& val);
  void splint(const Func& func, int k, double b, double step, double& val);
  void splint(const Func& func, int k, double b, double step, double& val,
              double& grad);

  // force calculation
  void doShift();  // make the rho / emf good
  void run();
  void run(int argc, char* argv[]);
  void calErr();
  void updaterho(vector<double>& vv);
  void updaterhoMEAM(vector<double>& vv);
  double errFunct(const vector<double>& x);
  double errFunctGrad(const vector<double>& x, vector<double>& g);
  double forceEAM(vector<Func>& ffs, int tag);
  double forceEAM(const vector<double>& vv);
  double forceADP(const vector<double>& vv, int tag);
  double forceMEAM(const vector<double>& vv);
  double forceEAM(const arma::mat& vv);
  double forceMEAM(const arma::mat& vv);
  void forceMEAM(Config& cc);
  void forceEAM(Config& cc);
  void stressMEAM(Config& cc);

  // haven't use it yet (need debug)
  double forceMEAMtb(const vector<double>& vv, int tag);
  void changePHI();
  void changeRHO();
  void changeEMF();
  void changeMEAMF();
  void changeMEAMG();

  // optimization
  void simAnneal();
  void simAnnealEAM();
  void randomize(vector<double>& vv, const int n, const vector<double>& v);
  int rescaleEMF(vector<double>& vv);
  int rescaleEMF();
  int rescaleRHO(vector<double>& vv);
  void shiftRHO(vector<double>& vv);
  void shiftEMF(double shift);
  void nloptGlobal();

  // diff evo
  void initPop();
  void diffEvo();

  // cmaes
  void loopcmaes();
  void cntcmaes();
  double cmaes(arma::mat& iterate);
  double testFunc(arma::mat& coordinates);

  // random
  double randNormal();
  double randUniform();
  double randUniform(const double min, const double max);

  // increase nodes
  void upgrade(int id);
  void increAnneal();
  void recordStage(int cnt);

  // inputs
  void readConfig();
  void readPot();
  void readParam();
  void readLmpMEAM();

  // outputs
  void writePot();
  void writePot(const vector<double>& vv);
  void writeLMPS();
  void writeLMPS(const vector<double>& vv);
  void writeMEAM();
  void writePOSCAR(const Config& cc, string fnm = "POSCAR.vasp");

  // utils
  void outMkdir(string mdir);
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

  void calLat(string key);
  void calPV();
  void calElas();
  void calSurf();

  // MPI utilitis
  void pfMPIinit(int argc, char* argv[]);
  void pfMPIfinalize();
  int createMPIdataTypes();

  Config addvolm(const Config&, const double& dl);
  Config addotho(const Config&, const double& dl);
  Config addmono(const Config&, const double& dl);
  Config addstrain(Config cc, const vector<vector<double>>& str);

  // analysize
  void deleteAtoms();
  void cutoffNeighs();
  void loopBwth();
  void forceDis();

  // debug
  void testSpline();
  friend class pfLMPdrv;
  friend class pfOptimizer;
};

class pfUtil {
 public:
  void split(const string& s, const char* delim, vector<string>& v);
  friend class pfHome;
};

#endif  // pfHome_H_
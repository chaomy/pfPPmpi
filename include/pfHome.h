#ifndef pfHome_H_
#define pfHome_H_

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
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

using std::string;
using std::unordered_map;
using std::vector;
namespace mpi = boost::mpi;
typedef enum { FCC, BCC, HCP, DAM, DIA, B1, C11, L12, B2 } lattice_t;

class pfOptimizer;
class pfHome {
 public:
  mpi::environment env;
  mpi::communicator cmm;
  mpi::communicator cmmlm;

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
  unordered_map<string, vector<double>> meamparms;
  vector<string> elems;     //  element name
  vector<lattice_t> lattp;  //  FCC, BCC, HCP, DAM, DIA, B1, C11, L12, B2
  vector<int> cnn1;         //  num near neighbors
  vector<double> t0;        //  1
  vector<int> ielement;     //  1
  vector<int> ibar;         //  3
  vector<double> rozero;
  vector<string> latticemp;
  vector<double> atwt;
  vector<double> alat;

  /* tests */
  Config ubcc;  // primitive bcc
  Config cbcc;  // conventional bcc
  unordered_map<string, vector<Config>> mpcf;
  unordered_map<string, vector<double>> mpvc;

  unordered_map<string, double> targs;
  unordered_map<string, double> exprs;
  unordered_map<string, double> weigh;
  unordered_map<string, double> error;

  /* map functions */
  unordered_map<string, void (pfHome::*)(Config&)> calfrc;
  unordered_map<string, double (pfHome::*)(const arma::mat& vv, int tg)> calobj;
  unordered_map<string, void (pfHome::*)()> write;

  vector<Config> configs;
  vector<Func> funcs;
  vector<Func> fprec;

  vector<int> startps;
  vector<int> endps;
  class pfLMPdrv;
  pfLMPdrv* lmpdrv;
  pfOptimizer* optdrv;
  Melem mele;

  vector<double> mfrc;
  vector<double> hil;  // 5 + 1
  vector<double> lol;  // 5 + 1
  vector<double> recorderr;
  vector<double> ini;
  vector<double> hib;
  vector<double> lob;
  vector<double> deb;  // hib - lob
  vector<int> gradRight;

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

  // void setupMEAMC();  // setup meam params
  void wrapAtomPos(Config& cc);
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
  void calErr(int tm);  // for debugging
  void updaterho(vector<double>& vv);
  void updaterhoMEAM(vector<double>& vv);
  double errFunct(const vector<double>& x);
  double errFunctGrad(const vector<double>& x, vector<double>& g);
  double forceEAM(vector<Func>& ffs, int tag);
  double forceEAM(const vector<double>& vv);
  double forceADP(const vector<double>& vv, int tag);
  double forceMEAM(const vector<double>& vv);
  double forceEAM(const arma::mat& vv);
  double forceEAM(const arma::mat& vv, int tg);
  double forceMEAM(const arma::mat& vv);
  double forceMEAM(const arma::mat& vv, int tg);
  double forceMEAMC(const arma::mat& vv, int tg);
  void forceMEAMC(Config& cc);
  void forceMEAM(Config& cc);
  void forceEAM(Config& cc);
  void stressMEAM(Config& cc);

  // optimization
  arma::mat encodev(const arma::mat& vv);
  arma::mat encodev(const vector<double>& vv);
  arma::mat decodev(const arma::mat& vv);
  vector<double> decodestdv(const vector<double>& vv);
  vector<double> encodestdv(const vector<double>& vv);
  void simAnneal();
  void simAnnealSpline();
  void randomize(vector<double>& vv, const int n, const vector<double>& v);
  void randomizeSpline(vector<double>& vv, const int n,
                       const vector<double>& v);
  int rescaleEMF(vector<double>& vv);
  int rescaleEMF(arma::mat& vv);
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
  void readMEAMC();
  void readParam();
  void readLmpMEAM();

  // outputs
  void writePot();
  void writePot(const vector<double>& vv);
  void writePot(const string& s);
  void writeLMPS();
  void writeLMPS(const vector<double>& vv);
  void writeMEAM();
  void writePOSCAR(const Config& cc, string fnm = "POSCAR.vasp");
  void writeMEAMC();

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
  void bcdata();

  Config addvolm(const Config&, const double& dl);
  Config addotho(const Config&, const double& dl);
  Config addmono(const Config&, const double& dl);
  Config addstrain(Config cc, const vector<vector<double>>& str);

  // analysize
  void deleteAtoms();
  void cutoffNeighs();
  void loopBwth();
  void forceDis();

  // MEAM
 private:
  int nelt;
  double cutforce, cutforcesq;
  vector<vector<double>> Ec_meam, re_meam;
  vector<double> Omega_meam, Z_meam, A_meam;
  vector<vector<double>> alpha_meam;
  vector<double> rho0_meam;
  vector<vector<double>> delta_meam;
  vector<double> beta0_meam, beta1_meam, beta2_meam, beta3_meam;
  vector<double> t0_meam, t1_meam, t2_meam, t3_meam, rho_ref_meam;
  vector<int> ibar_meam, ielt_meam;
  vector<vector<lattice_t>> lattce_meam;
  vector<vector<int>> nn2_meam, zbl_meam, eltind;
  vector<vector<double>> attrac_meam, repuls_meam;
  vector<vector<vector<double>>> Cmin_meam, Cmax_meam;
  vector<vector<double>> ebound_meam;

  vector<vector<double>> phir, phirar, phirar1, phirar2, phirar3, phirar4,
      phirar5, phirar6;
  double rc_meam, delr_meam, gsmooth_factor;
  int augt1, ialloy, mix_ref_t, emb_lin_neg, bkgd_dyn, erose_form;
  int vind2D[3][3], vind3D[3][3][3];
  int v2D[6], v3D[10];

  int nr, nrar;
  double dr, rdrar;

 public:
  int errorflag;
  vector<double> rho, rho0, rho1, rho2, rho3, frhop;
  vector<double> gamma, dgamma1, dgamma2, dgamma3, arho2b;
  vector<vector<double>> arho1, arho2, arho3, arho3b, t_ave, tsq_ave;

 private:
  static double fcut(const double xi);
  static double dfcut(const double xi, double& dfc);
  static double dCfunc(const double rij2, const double rik2, const double rjk2);
  static void dCfunc2(const double rij2, const double rik2, const double rjk2,
                      double& dCikj1, double& dCikj2);

  double G_gam(const double gamma, const int ibar, int& errorflag) const;
  double dG_gam(const double gamma, const int ibar, double& dG) const;
  static double zbl(const double r, const int z1, const int z2);
  static double erose(const double r, const double re, const double alpha,
                      const double Ec, const double repuls, const double attrac,
                      const int form);

  static void get_shpfcn(const lattice_t latt, double (&s)[3]);
  static int get_Zij(const lattice_t latt);
  static int get_Zij2(const lattice_t latt, const double cmin,
                      const double cmax, double& a, double& S);

 protected:
  void meam_checkindex(int, int, int, int*, int*);
  void getscreen(Config& cc);
  void calc_rho1(Config& cc);
  void alloyparams();
  void compute_pair_meam();
  void compute_pair_meam(int debug);
  double phi_meam(double, int, int);
  void compute_reference_density();
  void get_tavref(double*, double*, double*, double*, double*, double*, double,
                  double, double, double, double, double, double, int, int,
                  lattice_t);
  void get_sijk(double, int, int, int, double*);
  void get_densref(double, int, int, double*, double*, double*, double*,
                   double*, double*, double*, double*);
  void interpolate_meam(int);
  double compute_phi(double, int, int);

  void meam_setup_global();
  void meam_setup_globalfixed();
  void meam_setup_global(const vector<double>& vv);
  void meam_setup_global(const arma::mat& vv);
  void meam_setup_param(int which, double value, int nindex,
                        int* index /*index(3)*/, int* errorflag);

  void meam_setup_done();
  void meam_dens_setup(Config& cc);
  void meam_dens_init(Config& cc);
  void meam_dens_final(Config& cc);
  void meam_force(Config& cc);

  // debug
  void testSpline();
  friend class pfOptimizer;
};

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline arma::mat pfHome::decodev(const arma::mat& vv) {
  arma::mat rs(nvars, 1);
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

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline vector<double> pfHome::decodestdv(const vector<double>& vv) {
  vector<double> rs(vv.size());
  for (int i = 0; i < nvars; i++)
    rs[i] = lob[i] + deb[i] * 0.5 * (1. - cos(PI * 0.1 * vv[i]));
  return rs;
}

// [a, b] -> [0, 10]
inline vector<double> pfHome::encodestdv(const vector<double>& vv) {
  vector<double> rs(vv.size());
  for (int i = 0; i < nvars; i++)
    rs[i] = 10 * acos(1. - 2. / deb[i] * (vv[i] - lob[i])) * INVPI;
  return rs;
}

class pfUtil {
 public:
  void split(const string& s, const char* delim, vector<string>& v);
  friend class pfHome;
};

/* cutoff function (16) */
inline double pfHome::fcut(const double xi) {
  double a;
  if (xi >= 1.0)
    return 1.0;
  else if (xi <= 0.0)
    return 0.0;
  else {
    a = 1.0 - xi;
    a *= a;
    a *= a;
    a = 1.0 - a;
    return a * a;
  }
}

/* cutoff function and its derivative (16) */
inline double pfHome::dfcut(const double xi, double& dfc) {
  double a, a3, a4, a1m4;
  if (xi >= 1.0) {
    dfc = 0.0;
    return 1.0;
  } else if (xi <= 0.0) {
    dfc = 0.0;
    return 0.0;
  } else {
    a = 1.0 - xi;
    a3 = a * a * a;
    a4 = a * a3;
    a1m4 = 1.0 - a4;

    dfc = 8 * a1m4 * a3;
    return a1m4 * a1m4;
  }
}

//-----------------------------------------------------------------------------
// Derivative of Cikj w.r.t. rij
//     Inputs: rij,rij2,rik2,rjk2
//
inline double pfHome::dCfunc(const double rij2, const double rik2,
                             const double rjk2) {
  double rij4, a, asq, b, denom;

  rij4 = rij2 * rij2;
  a = rik2 - rjk2;
  b = rik2 + rjk2;
  asq = a * a;
  denom = rij4 - asq;
  denom = denom * denom;
  return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom;
}

//-----------------------------------------------------------------------------
// Derivative of Cikj w.r.t. rik and rjk
//     Inputs: rij,rij2,rik2,rjk2
//
inline void pfHome::dCfunc2(const double rij2, const double rik2,
                            const double rjk2, double& dCikj1, double& dCikj2) {
  double rij4, rik4, rjk4, a, denom;

  rij4 = rij2 * rij2;
  rik4 = rik2 * rik2;
  rjk4 = rjk2 * rjk2;
  a = rik2 - rjk2;
  denom = rij4 - a * a;
  denom = denom * denom;
  dCikj1 = 4 * rij2 *
           (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom;
  dCikj2 = 4 * rij2 *
           (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom;
}

static inline bool iszero(const double f) { return fabs(f) < 1e-20; }

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

#endif  // pfHome_H_
#ifndef PF_MEAMC_H
#define PF_MEAMC_H

#include "pfForce.h"
#include "pfHome.h"

#define MXEL 5  // max number of elements to handle

typedef enum { FCC, BCC, HCP, DAM, DIA, B1, C11, L12, B2 } lattice_t;

// MEAM Analytic form
class pfHome::pfForce::pfMEAMC {
  /* calculations follow the routines in LAMMPS  */
  vector<lattice_t> lattp;  //  FCC, BCC, HCP, DAM, DIA, B1, C11, L12, B2
  vector<int> cnn1;         //  num near neighbors
  vector<double> t0;        //  1
  vector<int> ielement;     //  1
  vector<int> ibar;         //  3
  vector<double> rozero;
  vector<string> latticemp;
  vector<double> atwt;
  vector<double> alat;

 private:
  pfHome& hm;
  mpi::communicator& cmm;
  unordered_map<string, string>& sparams;
  unordered_map<string, vector<double>>& meamparms;
  vector<double>& ini;
  vector<string>& elems;
  vector<double>& lob;
  vector<double>& hib;  // hold high bound spline node values
  vector<double>& deb;      // hib - lob
  vector<Config>& configs;  // configurations
  unordered_map<string, double>& dparams;
  unordered_map<string, int>& iparams;
  unordered_map<string, double>& error;
  int& nconfs;
  int& locstt;
  int& locend;
  int& nvars;

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

  pfMEAMC(pfHome& x)
      : hm(x),
        cmm(x.cmm),
        sparams(x.sparams),
        meamparms(x.meamparms),
        ini(x.ini),
        elems(x.elems),
        lob(x.lob),
        hib(x.hib),
        deb(x.deb),
        configs(x.configs),
        dparams(x.dparams),
        iparams(x.iparams),
        error(x.error),
        nconfs(x.nconfs),
        locstt(x.locstt),
        locend(x.locend),
        nvars(x.nvars),
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
        Cmin_meam(MXEL,
                  vector<vector<double>>(MXEL, vector<double>(MXEL, 0.49))),
        Cmax_meam(MXEL,
                  vector<vector<double>>(MXEL, vector<double>(MXEL, 2.8))),
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
    latticemp = vector<string>(
        {"fcc", "bcc", "hcp", "dim", "dia", "b1", "c11", "l12", "b2"});
  }

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

  void bcdata();
  void readMEAMC();
  void readMEAMCcnt();
  void writeMEAMC();
  double forceMEAMC(const arma::mat& vv, int tg);
  void forceMEAMC(Config& cc);

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
};

/* cutoff function (16) */
inline double pfHome::pfForce::pfMEAMC::fcut(const double xi) {
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
inline double pfHome::pfForce::pfMEAMC::dfcut(const double xi, double& dfc) {
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
inline double pfHome::pfForce::pfMEAMC::dCfunc(const double rij2,
                                               const double rik2,
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
inline void pfHome::pfForce::pfMEAMC::dCfunc2(const double rij2,
                                              const double rik2,
                                              const double rjk2, double& dCikj1,
                                              double& dCikj2) {
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

#endif  // PF_MEAMC_H
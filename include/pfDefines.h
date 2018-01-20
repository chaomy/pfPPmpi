#ifndef pfDefines
#define pfDefines

#include <vector>
#include "spline.h"
using std::cout;
using std::endl;
using std::string;
using std::vector;

#define DIM 3
#define MAXLEN 1024
#define PWEIGHT 100
#define EWEIGH 500
#define FRCEPS 0.1
#define PFROOT 0
#define PFERROR 1
#define PFSUCCESS 0
#define EVA3_GPA 160.21766208

enum {
  PHI = 0,
  RHO = 1,
  EMF = 2,
  ADPU = 3,
  ADPW = 4,
  ADPC = 5,
  LMPPNTS = 10000
};

// nearest neighbors
// 1st 8 ; 2nd 6; 3rd 12; 4th 24; 5th 8
enum { N1 = 8, N2 = 14, N3 = 26, N4 = 50, N5 = 58 };
enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };
enum { MEAMF = 3, MEAMG = 4 };
typedef double vec3[3];
typedef double position[3];
typedef double force[3];
typedef double vec6[6];

inline void printVec(double v[3]) {
  for (int i = 0; i < 3; i++) cout << v[i] << " ";
  cout << endl;
}

inline double square11(double x) { return x * x; };
inline double square33(double v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
};

inline double innDot33(double a[3], double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

inline double relerr(double a, double b) { return fabs(a - b) / b; }

inline void scaleVec(double v[3], double s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
};

class Angle {
 private:
  double gcos;
  int slot;
  double step;
  double shift;
  double gval;
  double ggrad;
  double apart;
  friend class pfAtom;
  friend class pfHome;
};

class Neigh {
 private:
  int aid;
  double invr;
  double psum;
  double r;
  double dist[3];
  double dist2r[3];
  double dis2mat[6];  // xx yy zz xy yz zx
  double phi, phig;
  double rho, rhog;
  double uval, ugrad;
  double wval, wgrad;
  double fval, fgrad;
  vector<int> slots;
  vector<double> steps;
  vector<double> shifts;
  friend class pfAtom;
  friend class pfHome;
};

typedef struct {
  int id;
  position pst;
  position prl;
  force frc;
  double absfrc;
} AtomMPI;

class pfAtom {
 private:
  int id;
  int tp;
  int nneighs;
  int nneighsFull;
  position pst;
  position prl;

  force frc;
  force fitfrc;
  force fweigh;

  force phifrc;
  force rhofrc;
  force trifrc;

  /* ADP */
  double mu[3];
  double nu;
  double lambda[6];

  vector<Neigh> neighs;

  /* MEAM */
  vector<Neigh> neighsFull;
  vector<vector<Angle>> angMat;

  double absfrc;
  double gradF;
  double rho;
  double crho;
  double prho;

 public:
  pfAtom() : id(0) {
    fitfrc[0] = 0.0;
    fitfrc[1] = 0.0;
    fitfrc[2] = 0.0;
  };
  pfAtom(int n) : id(n) {
    fitfrc[0] = 0.0;
    fitfrc[1] = 0.0;
    fitfrc[2] = 0.0;
  };
  ~pfAtom(){};
  friend class Config;
  friend class pfHome;
};

class Config {
 private:
  int natoms;
  int cfgid;
  double weigh;
  double engy;
  double fitengy;
  double phiengy;
  double emfengy;
  double vol;
  double rhomx;
  double rhomi;
  int scale[3];
  vec3 bvx, tvx;
  vec3 bvy, tvy;
  vec3 bvz, tvz;
  vec6 strs;

  vector<pfAtom> atoms;
  vector<int> natomsv;
  vector<string> nelemsv;

 public:
  Config() { cfgid = 0, natoms = 0, weigh = 0.0, engy = 0.0; };
  Config(int n) : cfgid(n) { natoms = 0, weigh = 0.0, engy = 0.0; };
  ~Config(){};
  friend class pfHome;
};

class Func {
 private:
  int npts;
  int bnd;
  double rng;      // range ymax - ymin
  double step;     // equidist; has't initilaized yet
  double invstep;  // equidist
  vector<double> xx;
  vector<double> yy;
  vector<double> g1;
  vector<double> g2;
  tk::spline s;
  friend class pfHome;
};

#endif
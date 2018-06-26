#ifndef pfLMP_H_
#define pfLMP_H_

#include "atom.h"
#include "input.h"
#include "lammps.h"  // these are LAMMPS include files
#include "library.h"
#include "pfHome.h"

#define LMPROOT 0
#define EV_A2toJ_M2 16.0218
#define CLEN 10

using namespace LAMMPS_NS;

class pfHome::pfLMPdrv {
 private:
  // commonly used
  char ETOTAL[CLEN] = "etotal";
  char PRESS[CLEN] = "press";
  char LX[CLEN] = "lx";
  char LY[CLEN] = "ly";
  char LZ[CLEN] = "lz";
  char NATOM[CLEN] = "atoms";
  char CLEAR[CLEN] = "clear";
  char SXX[CLEN] = "Sxx";
  char SYY[CLEN] = "Syy";
  char SZZ[CLEN] = "Szz";
  char SXY[CLEN] = "Sxy";
  char SXZ[CLEN] = "Sxz";
  char SYZ[CLEN] = "Syz";

  int cnt = 0;

  MPI_Comm lmp_comm;
  int mrank, nprocs;
  LAMMPS* lmp = NULL;
  pfHome& hm;
  unordered_map<string, double>& dparams;
  unordered_map<string, string>& sparams;  // string parameters
  unordered_map<string, double>& targs;
  unordered_map<string, double>& exprs;
  unordered_map<string, double>& weigh;
  Ctrib& error;
  vector<string>& elems;
  mpi::communicator& cmmlm;
  Melem& mele;  // element data

  // itensile
  vector<vector<double>> bsm;
  unordered_map<string, vector<double>> iegyv;
  unordered_map<string, vector<double>> isxxv;

 public:
  vector<int> gsfpnts;
  unordered_map<string, int> label;
  unordered_map<string, vector<double>> lgsf;  // gsf
  unordered_map<string, vector<double>> lmpv;  // pv
  pfLMPdrv(int argc, char* argv[], pfHome& x);
  ~pfLMPdrv();

  void calPhy();
  void calPhyErr();

  void calLatticeBCC();
  void calLatticeHCP();
  void calLatticeFCC();

  void calElastic();
  void calSurface();
  void calSurfaceUrlx();
  void calGSF();
  void calGSFUrlx();
  void calPV();
  void calVac();

  void calIten();
  void calItenOptLin(const double& dlt, const string& tag);
  void calItenRun(const string& tag, double& v, vector<double>& stsv);

  void paraInit(int argc, char* argv[]);
};

#endif  // pfLMP_H_
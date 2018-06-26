/*
 * @Author: chaomy
 * @Date:   2017-11-13 15:58:23
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-26 15:44:15
 */

#include "pfLmpDrv.h"

using namespace LAMMPS_NS;

pfHome::pfLMPdrv::pfLMPdrv(int argc, char* argv[], pfHome& x)
    : hm(x),
      dparams(x.dparams),
      sparams(x.sparams),
      targs(x.targs),
      exprs(x.exprs),
      weigh(x.weigh),
      perr(x.perr),
      elems(x.elems),
      cmmlm(x.cmmlm),
      mele(x.mele),
      gsfpnts(vector<int>({0, 3, 4, 5, 6, 7})) {
  paraInit(argc, argv);
}
pfHome::pfLMPdrv::~pfLMPdrv() {
  if (lmp) delete lmp;
}

void pfHome::pfLMPdrv::calPhy() {
  if (label["bcc"]) calLatticeBCC();
  if (label["fcc"]) calLatticeFCC();
  if (label["hcp"]) calLatticeHCP();
  if (label["cij"]) calElastic();
  if (label["suf"]) calSurface();
  if (label["gsf"]) calGSF();
  if (label["pv"]) calPV();
  // if (label["vac"]) calVac();
  // if (label["iten"]) calIten();
  calPhyErr();
  // remove("no");
  // remove("log.lammps");

  // exprs["bcc2hcp"] = exprs["ehcp"] - exprs["ebcc"];
  // exprs["bcc2fcc"] = exprs["efcc"] - exprs["ebcc"];
  // vector<string> aa({"lat", "bcc2fcc", "bcc2hcp", "c11", "c12", "c44",
  // "suf110",
  //                    "suf100", "suf111"});
  // for (string ee : aa)
  //   cout << ee << " " << exprs[ee] << " " << targs[ee] << endl;
}

void pfHome::pfLMPdrv::calPhyErr() {
  exprs["bcc2hcp"] = exprs["ehcp"] - exprs["ebcc"];
  exprs["bcc2fcc"] = exprs["efcc"] - exprs["ebcc"];
  vector<string> aa({"lat", "c11", "c12", "c44", "suf110", "suf100", "suf111",
                     "bcc2fcc", "bcc2hcp"});
  for (string ee : aa)
    cout << ee << " " << exprs[ee] << " " << targs[ee] << endl;
  for (string ee : aa) perr[ee] = 100 * relerr(exprs[ee], targs[ee]);
  for (string ee : aa) cout << "perr " << ee << " " << perr[ee] << endl;

  if (label["gsf"]) {
    perr["gsf110"] = perr["gsf211"] = 0.0;
    for (double ee : lgsf["111e110"]) perr["gsf110"] += ee;
    perr["gsf110"] /= lgsf["111e110"].size();

    for (double ee : lgsf["111e211"]) perr["gsf211"] += ee;
    perr["gsf211"] /= lgsf["111e211"].size();

    perr["gsf110"] *= 100;
    perr["gsf211"] *= 100;
  }

  // usf
  // exprs["usf211"] = exprs["usf110"] = 0.0;
  // for (double ee : lgsf["111z110"]) exprs["usf110"] = fmax(exprs["usf110"],
  // ee); for (double ee : lgsf["111z211"]) exprs["usf211"] =
  // fmax(exprs["usf211"], ee);

  // pv
  if (label["pv"]) {
    perr["pv"] = 0.0;
    for (double ee : lmpv["Nbe"]) perr["pv"] += ee;
    perr["pv"] /= lmpv["Nbe"].size();
    perr["pv"] *= 100;
  }

  if (label["gsf"] && label["pv"]) {  // output
    string fmt(
        "%3d %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f "
        "%06.2f\n");
    FILE* fid = fopen("perr.log", "a");
    fprintf(fid, fmt.c_str(), cnt++, exprs["lat"], exprs["c11"], exprs["c12"],
            exprs["c44"], exprs["suf110"], exprs["suf100"], exprs["suf111"],
            exprs["usf110"], exprs["usf211"], lmpv["Nbp"].back());
    fprintf(fid, fmt.c_str(), 0, perr["lat"], perr["c11"], perr["c12"],
            perr["c44"], perr["suf110"], perr["suf100"], perr["suf111"],
            perr["gsf110"], perr["gsf211"], perr["pv"]);
    fclose(fid);
  }
}

void pfHome::pfLMPdrv::paraInit(int argc, char* argv[]) {
  label["bcc"] = 1;
  label["fcc"] = 1;
  label["hcp"] = 1;
  label["cij"] = 1;
  label["suf"] = 1;
  label["gsf"] = 0;
  label["pv"] = 0;
  // label["vac"] = 1;
  // label["iten"] = 1;

  // mpi version without boost
  // MPI_Init(&argc, &argv);
  // MPI_Comm_rank(MPI_COMM_WORLD, &mrank);
  // MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  char** lmparg = new char*[3];
  lmparg[0] = NULL;  // placeholder for program name
  lmparg[1] = (char*)"-screen";
  lmparg[2] = (char*)"no";

  // lmp = new LAMMPS(3, lmparg, MPI_COMM_WORLD);
  lmp = new LAMMPS(3, lmparg, cmmlm);
  delete[] lmparg;
}

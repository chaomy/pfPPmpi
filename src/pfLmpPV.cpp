/*
 * @Author: chaomy
 * @Date:   2017-11-23 07:29:33
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-16 00:18:50
 */

#include "pfLmpDrv.h"

using std::pow;

void pfHome::pfLMPdrv::calPV() {
  string ptg(sparams["elem"] + "p");
  string vtg(sparams["elem"] + "v");
  string etg(sparams["elem"] + "e");

  const vector<double>& pp = mele.pvm[ptg];
  const vector<double>& vv = mele.pvm[vtg];

  for (int it = 0; it < pp.size(); it++) {
    char cmds[100][MAXLEN];
    int i = 0;
    //  --------------------- INITIALIZAITION ---------------------
    sprintf(cmds[i++], "clear");
    // sprintf(cmds[i++], "print logfile screen no");
    sprintf(cmds[i++], "units  metal");
    sprintf(cmds[i++], "dimension  3");
    sprintf(cmds[i++], "boundary p p p");
    sprintf(cmds[i++], "atom_style atomic");
    sprintf(cmds[i++], "variable  a equal  %.7f",
            pow(vv[it], 1. / 3.) * exprs["lat"]);

    // --------------------- ATOM DEFINITION ---------------------
    sprintf(cmds[i++], "lattice  bcc  ${a}");
    sprintf(cmds[i++], "region whole block 0 ${a} 0 ${a} 0 ${a}  units  box");
    sprintf(cmds[i++], "create_box  1  whole");
    sprintf(
        cmds[i++],
        "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 0  orient  z 0 0 1");
    sprintf(cmds[i++], "create_atoms 1 region whole");

    // --------------------- FORCE FIELDS ---------------------
    sprintf(cmds[i++], "pair_style  %s", sparams["pairstyle"].c_str());
    if (!sparams["ptype"].compare("MEAMC"))
      sprintf(cmds[i++], "pair_coeff  *  *  %s %s %s %s",
              sparams["meamlib"].c_str(), elems[0].c_str(),
              sparams["meampar"].c_str(), elems[0].c_str());
    else
      sprintf(cmds[i++], "pair_coeff * * %s %s", sparams["lmpfile"].c_str(),
              sparams["elem"].c_str());

    // sprintf(cmds[i++], "mass  *  %f", dparams["mass"]);
    sprintf(cmds[i++], "neighbor 1.0 bin");
    sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

    // ---------------------- THERMO  --------------------------
    sprintf(cmds[i++], "thermo 10000");
    sprintf(cmds[i++],
            "thermo_style  custom  lx  ly  lz step  etotal temp press");
    // ---------------------- RELAX   --------------------------
    sprintf(cmds[i++], "min_style  cg");
    sprintf(cmds[i++], "minimize  1e-12  1e-12  100000  100000");

    for (int mm = 0; mm < i; mm++) lammps_command(lmp, cmds[mm]);

    lmpv[ptg].push_back(1e-4 * lammps_get_thermo(lmp, PRESS));

    double ee;
    if (pp[it] != 0)
      ee = abs(lmpv[ptg][it] - pp[it]) / pp[it];
    else
      ee = abs(lmpv[ptg][it] - pp[it]) / 1.0;

    lmpv[etg].push_back(ee);

    lammps_command(lmp, CLEAR);
    printf("%d %f %f %f %f\n", it, vv[it], pp[it], lmpv[ptg][it],
           lmpv[etg][it]);
  }
}

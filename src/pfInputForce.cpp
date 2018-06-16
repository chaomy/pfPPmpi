/*
 * @Author: chaomy
 * @Date:   2018-01-20 16:53:38
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-16 16:01:41
 */

#include "pfIO.h"

void pfHome::pfIO::readConfig() { /* read atomic force file */
  double eb = -3.09477080931;     // Energy per atom shift -3.09477080931
  // double eperf = -6.9950724562545;
  double eperf = -7.01;
  configs.clear();  // clear
  char tmp[MAXLEN];

  ifstream ifs(sparams["cnffile"].c_str(), std::ifstream::in);

  int cnt = 0;
  string buff;
  vector<string> segs;

  Config cnf;
  while (getline(ifs, buff)) {
    segs.clear();
    cnf.atoms.clear();
    split(buff, " ", segs);
    if (!segs[0].compare("#N")) {
      sscanf(buff.c_str(), "%s %d %s", tmp, &cnf.natoms, tmp);
    } else if (!segs[0].compare("#X")) {
      sscanf(buff.c_str(), "%s %lf %lf %lf", tmp, &cnf.bvx[0], &cnf.bvx[1],
             &cnf.bvx[2]);
    } else if (segs[0] == "#Y") {
      sscanf(buff.c_str(), "%s %lf %lf %lf", tmp, &cnf.bvy[0], &cnf.bvy[1],
             &cnf.bvy[2]);
    } else if (segs[0] == "#Z") {
      sscanf(buff.c_str(), "%s %lf %lf %lf", tmp, &cnf.bvz[0], &cnf.bvz[1],
             &cnf.bvz[2]);
    } else if (!segs[0].compare("#W")) {
      sscanf(buff.c_str(), "%s %lf", tmp, &cnf.weigh);
    } else if (!segs[0].compare("#E")) {
      sscanf(buff.c_str(), "%s %lf", tmp, &cnf.engy);
      cnf.engy -= eb;
    } else if (!segs[0].compare("#S")) {
      sscanf(buff.c_str(), "%s %lf %lf %lf %lf %lf %lf", tmp, &cnf.strs[0],
             &cnf.strs[1], &cnf.strs[2], &cnf.strs[3], &cnf.strs[4],
             &cnf.strs[5]);
    } else if (!segs[0].compare("#F")) {
      for (int i = 0; i < cnf.natoms; i++) {
        getline(ifs, buff);
        pfAtom atom(i);
        sscanf(buff.c_str(), "%s %lf %lf %lf %lf %lf %lf", tmp, &atom.pst[0],
               &atom.pst[1], &atom.pst[2], &atom.frc[0], &atom.frc[1],
               &atom.frc[2]);
        atom.absfrc = sqrt(square33(atom.frc));
        for (int ii : {0, 1, 2}) {
          atom.fweigh[ii] = 1. / (fabs(atom.frc[ii]) + dparams["frceps"]);
          // atom.fweigh[ii] =
          //     1e4 * exp(-dparams["bwidth"] * atom.frc[ii] * atom.frc[ii]);
          // mfrc[ii] += fabs(atom.frc[ii]);
        }
        cnf.atoms.push_back(atom);
      }
      tln += cnf.natoms;
      cnf.cfgid = cnt++;
      cnf.weigh =
          cnf.weigh * std::exp(-fabs(cnf.engy - eperf) / dparams["bwidth"]);
      configs.push_back(cnf);
    }  // #F
  }    // while

  ifs.close();
  cout << "finish reading " << configs.size() << " configs" << endl;
  nconfs = configs.size();
}
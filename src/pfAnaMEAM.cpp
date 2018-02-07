/*
 * @Author: chaomy
 * @Date:   2018-01-29 22:10:28
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-06 19:02:18
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

double pfHome::forceMEAMC(const arma::mat& vv, int tg) {
  while (true) {
    broadcast(cmm, tg, PFROOT);
    if (tg == EXT) break;

    for (int i = 0; i < nvars; i++) ini[i] = vv[i];
    broadcast(cmm, ini, PFROOT);

    meam_setup_global(ini);
    meam_setup_done();
    double efrc = 0.0;
    for (int i = locstt; i < locend; i++) {
      Config& cc = configs[i];
      forceMEAMC(cc);
      for (pfAtom& atm : cc.atoms) {
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] = -atm.frc[it];
          efrc += square11(atm.fitfrc[it] * atm.fweigh[it]);
        }
      }
      cc.fitengy /= cc.natoms;
      efrc += square11(cc.fitengy - cc.engy);
    }
    reduce(cmm, efrc, error["frc"], std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  return error["frc"];
}

void pfHome::forceMEAMC(Config& cc) {
  cc.fitengy = 0.0;
  for (pfAtom& atm : cc.atoms) {
    atm.eng = 0.0;
    for (int ii : {0, 1, 2}) atm.fitfrc[ii] = 0.;
    for (int ii : {0, 1, 2, 3, 4, 5}) atm.sts[ii] = 0.;
  }
  meam_dens_setup(cc);
  meam_dens_init(cc);
  meam_dens_final(cc);
  meam_force(cc);
}

// For Debugging
// if (cmm.rank() == 1) {
//   FILE* fid = fopen("pf.txt", "w");
//   for (int i = 0; i < cc.natoms; i++) {
//     pfAtom& atm = cc.atoms[i];
//     fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f %.4f \n", atm.pst[0],
//             atm.pst[1], atm.pst[2], atm.fitfrc[0], atm.fitfrc[1],
//             atm.fitfrc[2], atm.eng);
//   }
//   fclose(fid);
// }

// FILE* fid = fopen("pf.txt", "w");
// for (int i = 0; i < 200; i++) {
//   vector<vector<double>> vc;
//   pfAtom& atm = cc.atoms[i];
//   for (int j : atm.neighidxHalf) {
//     Neigh& ngbj = atm.neighsFull[j];
//     vc.push_back(vector<double>(
//         {atm.pst[0] + ngbj.dist[0], atm.pst[1] + ngbj.dist[1],
//          atm.pst[2] + ngbj.dist[2], ngbj.dscrfcn, ngbj.scrfcn,
//          ngbj.fcpair}));
//   }
//   sort(vc.begin(), vc.end());
//   for (auto ee : vc)
//     fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f\n", ee[0], ee[1], ee[2],
//             ee[3], ee[4], ee[5]);

//   vector<double> tm({atm.pst[0], atm.pst[1], atm.pst[2], atm.fitfrc[0],
//                      atm.fitfrc[1], atm.fitfrc[2], atm.eng});
//   fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", tm[0], tm[1], tm[2],
//           tm[3], tm[4], tm[5], tm[6]);
// }

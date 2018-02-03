/*
 * @Author: chaomy
 * @Date:   2018-01-29 22:10:28
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-02 23:45:11
 */

#include "pfHome.h"

double pfHome::forceMEAMC() {
  Config& cc = configs[0];
  for (pfAtom& atm : cc.atoms) {
    for (int ii : {0, 1, 2}) atm.fitfrc[ii] = 0.;
    for (int ii : {0, 1, 2, 3, 4, 5}) atm.sts[ii] = 0.;
  }

  meam_dens_setup(cc);
  meam_dens_init(cc);
  meam_dens_final(cc);
  meam_force(cc);

  FILE* fid = fopen("pf.txt", "w");
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

  for (int i = 0; i < cc.natoms; i++) {
    pfAtom& atm = cc.atoms[i];
    fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f %.4f \n", atm.pst[0],
            atm.pst[1], atm.pst[2], atm.fitfrc[0], atm.fitfrc[1], atm.fitfrc[2],
            atm.eng);
  }

  fclose(fid);
  return 0;
}

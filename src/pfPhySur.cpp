/*
 * @Author: chaomy
 * @Date:   2017-12-19 16:23:55
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-16 17:00:30
 */

#include "pfConf.h"
#include "pfForce.h"
#include "pfPhy.h"

void pfHome::pfPhy::calSurf(pfForce& fcdrv, pfConf& cfdrv) {
  double la = exprs["lat"];
  printf("la = %f\n", la);
  Config s1 = cfdrv.buildsur100(la, "sur");
  (fcdrv.*calfrc[sparams["ptype"]])(s1);
  printf("%f\n",
         (s1.fitengy - cfdrv.ubcc.fitengy) * 0.5 * s1.natoms / (la * la));
}

Config pfHome::pfConf::buildsur100(const double& lat,
                                   const string& tag) { /* 100 */
  Config cc;
  int h = 24;
  int l = (!tag.compare("sur")) ? h - 10 : h;
  double lz = h * lat;

  cc.bvx[X] = lat, cc.bvx[Y] = 0.0, cc.bvx[Z] = 0.0;
  cc.bvy[X] = 0.0, cc.bvy[Y] = lat, cc.bvy[Z] = 0.0;
  cc.bvz[X] = 0.0, cc.bvz[Y] = 0.0, cc.bvz[Z] = lz;

  int cn = 0;
  for (int iz = 0; iz < l; iz++) {
    pfAtom atm1(cn++);
    atm1.prl[X] = atm1.prl[Y] = 0.0;
    atm1.prl[Z] = double(iz) / double(h);

    atm1.pst[X] = atm1.prl[X] * lat;
    atm1.pst[Y] = atm1.prl[Y] * lat;
    atm1.pst[Z] = atm1.prl[Z] * lz;
    cc.atoms.push_back(atm1);

    pfAtom atm2(cn++);
    atm2.prl[X] = atm2.prl[Y] = 0.5;
    atm2.prl[Z] = double(iz + 0.5) / double(h);

    atm2.pst[X] = atm2.prl[X] * lat;
    atm2.pst[Y] = atm2.prl[Y] * lat;
    atm2.pst[Z] = atm2.prl[Z] * lz;
    cc.atoms.push_back(atm2);
  }
  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}

Config pfHome::pfConf::buildsur110(const double& lat, const string& tag) {
  Config cc;
  return cc;
}

Config pfHome::pfConf::buildsur211(const double& lat, const string& tag) {
  Config cc;
  return cc;
}
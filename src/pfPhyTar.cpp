/*
 * @Author: chaomy
 * @Date:   2017-11-14 14:24:20
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-26 16:07:26
 */

#include "pfHome.h"

void pfHome::initTargs() {
  perr["tol"] = 0.0;

  // lattice
  // targs["lat"] = 3.308;
  targs["lat"] = 3.313;
  targs["abcc"] = 3.313;
  targs["afcc"] = 4.220;
  targs["ahcp"] = 2.90;
  targs["chcp"] = 5.27;

  // engy
  targs["ebcc"] = -10.0898;
  targs["efcc"] = -9.7707;
  targs["ehcp"] = -9.7952;

  targs["bcc2hcp"] = 0.295;
  targs["bcc2fcc"] = 0.320;

  // elastic
  targs["c11"] = 250.;
  targs["c12"] = 133.;
  targs["c44"] = 30.;

  // surf
  targs["suf100"] = 2.35;
  targs["suf110"] = 2.10;
  targs["suf111"] = 2.40;

  pwgh["lat"] = 1e2;
  pwgh["c11"] = 10;
  pwgh["c12"] = 10;
  pwgh["c44"] = 10;
  pwgh["suf"] = 1.;
}
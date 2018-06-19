/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:11:45
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-18 18:13:38
 */

#include "pfHome.h"
#include "pfMEAMC.h"

void pfHome::bcdata() {
  /* params */
  broadcast(cmm, dparams, PFROOT);
  broadcast(cmm, iparams, PFROOT);
  broadcast(cmm, sparams, PFROOT);

  /* configurations */
  broadcast(cmm, nconfs, PFROOT);
  broadcast(cmm, configs, PFROOT);

  /* funcs */
  broadcast(cmm, nfuncs, PFROOT);
  broadcast(cmm, ini, PFROOT);
  broadcast(cmm, optidx, PFROOT);
  broadcast(cmm, smthidx, PFROOT);
  nvars = ini.size();
  broadcast(cmm, lob, PFROOT);
  broadcast(cmm, hib, PFROOT);
  broadcast(cmm, deb, PFROOT);
  broadcast(cmm, startps, PFROOT);
  broadcast(cmm, endps, PFROOT);

  if (cmm.rank() != PFROOT) {
    for (int i = 0; i < nfuncs; i++) {
      Func tmp;
      funcs.push_back(tmp);
    }
  }

  broadcast(cmm, nvars, PFROOT);
  for (int i = 0; i < nfuncs; i++) {
    broadcast(cmm, funcs[i].bl, PFROOT);
    broadcast(cmm, funcs[i].br, PFROOT);
    broadcast(cmm, funcs[i].npts, PFROOT);
    broadcast(cmm, funcs[i].xx, PFROOT);
    broadcast(cmm, funcs[i].yy, PFROOT);
    broadcast(cmm, funcs[i].g1, PFROOT);
    broadcast(cmm, funcs[i].g2, PFROOT);
    broadcast(cmm, funcs[i].rlxid, PFROOT);
  }

  // assign cutoffs
  ricut = funcs[PHI].xx.front();
  rocut = funcs[PHI].xx.back();
  rhcut = funcs[RHO].xx.back();

  vector<tk::spline::bd_type> bdmp(
      {tk::spline::first_deriv, tk::spline::second_deriv});
  vector<bool> odmp({true, true});  // set 1 and 2 both to be enforce linear

  for (Func& ff : funcs) {
    double al = (ff.bl == 1) ? ff.g1.front() : ff.g2.front();  // left bound
    double bl = (ff.br == 1) ? ff.g1.back() : ff.g2.back();    // right bound
    ff.s.set_boundary(bdmp[ff.bl - 1], al, bdmp[ff.br - 1], bl,
                      odmp[ff.br - 1]);
  }
}

void pfHome::pfForce::pfMEAMC::bcdata() {
  broadcast(cmm, dparams, PFROOT);
  broadcast(cmm, iparams, PFROOT);
  broadcast(cmm, sparams, PFROOT);

  /* configurations */
  broadcast(cmm, nconfs, PFROOT);
  broadcast(cmm, configs, PFROOT);

  broadcast(cmm, elems, PFROOT);
  broadcast(cmm, cnn1, PFROOT);
  broadcast(cmm, ielement, PFROOT);
  broadcast(cmm, atwt, PFROOT);
  broadcast(cmm, alat, PFROOT);
  broadcast(cmm, rozero, PFROOT);
  broadcast(cmm, ibar, PFROOT);
  broadcast(cmm, ini, PFROOT);
  broadcast(cmm, t0, PFROOT);
  nvars = ini.size();
  meam_setup_globalfixed();
}

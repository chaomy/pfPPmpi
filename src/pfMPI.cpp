/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:11:45
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-07 21:14:47
 */

#include "pfHome.h"

using std::cout;
using std::endl;

void pfHome::bcdata() {
  /* params */
  broadcast(cmm, dparams, PFROOT);
  broadcast(cmm, iparams, PFROOT);
  broadcast(cmm, sparams, PFROOT);

  /* configurations */
  broadcast(cmm, nconfs, PFROOT);
  broadcast(cmm, configs, PFROOT);

  /* funcstions */
  if (!sparams["ptype"].compare("MEAMC")) {
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
  } else {  // use spline
    broadcast(cmm, nfuncs, PFROOT);
    broadcast(cmm, ini, PFROOT);
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
      broadcast(cmm, funcs[i].npts, PFROOT);
      broadcast(cmm, funcs[i].xx, PFROOT);
      broadcast(cmm, funcs[i].yy, PFROOT);
      broadcast(cmm, funcs[i].g1, PFROOT);
      broadcast(cmm, funcs[i].g2, PFROOT);
    }

    // assign cutoffs
    ricut = funcs[PHI].xx.front();
    rocut = funcs[PHI].xx.back();
    rhcut = funcs[RHO].xx.back();

    // cutforce = rocut;
    // cutforcesq = rocut * rocut;

    if (!sparams["ptype"].compare("MEAMS")) {  // boundary
      funcs[MEAMF].s.set_boundary(tk::spline::second_deriv, 0.0,
                                  tk::spline::first_deriv, 0.0, true);
      // funcs[MEAMG].s.set_boundary(tk::spline::second_deriv, 0.0,
      //                             tk::spline::second_deriv, 0.0, true);
      // funcs[MEAMF].s.set_boundary(
      //     tk::spline::first_deriv, funcs[MEAMF].g1.front(),
      //     tk::spline::first_deriv, funcs[MEAMF].g1.back(), true);
      funcs[MEAMG].s.set_boundary(
          tk::spline::first_deriv, funcs[MEAMG].g1.front(),
          tk::spline::first_deriv, funcs[MEAMG].g1.back(), true);
    }
    // funcs[PHI].s.set_boundary(tk::spline::second_deriv, 0.0,
    //                           tk::spline::first_deriv, 0.0, true);
    // funcs[RHO].s.set_boundary(tk::spline::second_deriv, 0.0,
    //                           tk::spline::first_deriv, 0.0, true);
    // funcs[EMF].s.set_boundary(tk::spline::second_deriv, 0.0,
    //                           tk::spline::second_deriv, 0.0, true);

    funcs[PHI].s.set_boundary(tk::spline::first_deriv, funcs[PHI].g1.front(),
                              tk::spline::first_deriv, funcs[PHI].g1.back(),
                              true);
    funcs[RHO].s.set_boundary(tk::spline::first_deriv, funcs[RHO].g1.front(),
                              tk::spline::first_deriv, funcs[RHO].g1.back(),
                              true);
    funcs[EMF].s.set_boundary(tk::spline::first_deriv, funcs[EMF].g1.front(),
                              tk::spline::first_deriv, funcs[EMF].g1.back(),
                              true);
  }
}

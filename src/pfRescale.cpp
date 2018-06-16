/*
 * @Author: chaomy
 * @Date:   2017-10-30 21:34:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 22:32:03
 */

#define NNEIGH 26
#include "pfForce.h"

// unfinished
void pfHome::shiftRHO(vector<double>& vv) {
  double shift = (0 - ominrho) / NNEIGH;
  for (unsigned int i = 0; i < funcs[RHO].yy.size(); i++)
    vv[startps[RHO] + i] += shift;
}

void pfHome::shiftEMF(double shift) {
  for (unsigned int i = 0; i < funcs[EMF].yy.size(); i++)
    funcs[EMF].xx[i] += shift;
}

int pfHome::rescaleEMF(arma::mat& in) {
  Func& ff = funcs[EMF];
  int npts = ff.npts;
  double delt = ff.step;
  if (fabs(ominrho - ff.xx.front()) > 2. * delt ||
      fabs(omaxrho - ff.xx.back()) > 2. * delt) {
    arma::mat iterate = decodev(in);  // -> [a, b]
    double ndelt = (omaxrho - ominrho) / (npts - 1);
    for (int i = 0; i < npts; i++) {
      double tt = ff.s(ff.xx[i] = ominrho + i * ndelt);
      int ii = startps[EMF] + i;
      if (tt < lob[ii])
        iterate[ii] = lob[ii] + 0.01 * arma::randu();
      else if (tt > lob[ii])
        iterate[ii] = hib[ii] - 0.01 * arma::randu();
      else
        iterate[ii] = tt;
    }
    ff.step = ndelt;
    in = encodev(iterate);
    return 1;
  } else
    return 0;
}

// REWRITE THE FUNCTION BEFORE USE
int pfHome::rescaleEMF(vector<double>& vv,
                       pfForce& fcdrv) {  // check if rescale
  Func& ff = funcs[EMF];
  int npts = ff.npts;
  double delt = ff.step;

  if (fabs(ominrho - ff.xx.front()) > delt ||
      fabs(omaxrho - ff.xx.back()) > delt) {
    for (int i = 0; i < npts; i++) ff.yy[i] = vv[startps[EMF] + i];

    ff.g1.front() = vv[nvars - 4];
    ff.g1.back() = vv[nvars - 2];

    double ndelt = (omaxrho - ominrho) / (npts - 1);

    /* update y values (to vv) */
    for (int i = 0; i < npts; i++)
      fcdrv.splint(ff, ominrho + i * ndelt, vv[startps[EMF] + i]);

    /* update xx values */
    for (int i = 0; i < npts; i++) ff.xx[i] = ominrho + i * ndelt;
    ff.step = ndelt;
    return (1);
  } else
    return (0);
}

int pfHome::rescaleRHO(vector<double>& vv) {
  printf("max %f hilim %f min %f lolim %f oaverho %f\n", omaxrho, hirho,
         ominrho, lorho, oaverho);
  double aa = (omaxrho - ominrho) / (hirho - lorho);
  if (aa > 1.2 || aa < 1.00) {  // rescale + shift
    double coeff = 1.01 / (aa);
    // rescale
    for (unsigned int i = 0; i < funcs[RHO].yy.size() + 2; i++)
      vv[startps[RHO] + i] *= coeff;
    // shift
    double shift = -coeff * (omaxrho + ominrho - hirho - lorho) / (2 * NNEIGH);
    for (unsigned int i = 0; i < funcs[RHO].yy.size(); i++)
      vv[startps[RHO] + i] += shift;
    printf("a = %f ; rescale %f ; shift %f \n", aa, coeff, shift);

  } else if ((omaxrho - hirho) * (ominrho - lorho) > 0) {
    double shift = -(omaxrho + ominrho - hirho - lorho) / (2 * NNEIGH);
    for (unsigned int i = 0; i < funcs[RHO].yy.size(); i++)
      vv[startps[RHO] + i] += shift;
    printf("a = %f ; shift %f \n", aa, shift);
  } else
    return (0);
  return (1);
}
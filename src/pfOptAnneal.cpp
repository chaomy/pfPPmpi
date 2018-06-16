/*
 * @Author: chaomy
 * @Date:   2017-10-23 20:10:54
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 22:44:07
 */

#include "pfForce.h"
#include "pfLmpDrv.h"

// #define EPS 0.1
#define EPS 10
#define TEMPVAR 0.85
#define STEPVAR 2.0
#define KMAX 100  // max # of runs  500 or 1000
#define INVSQRT2PI 0.39894228040143267794
#define GAUSS(a) (INVSQRT2PI * (exp(-((a) * (a)) / 2.0)))
#define NEPS 4
#define NSTEP 3  // 20
#define NTEMP (3 * 40)
#define RESFREQ 10

void pfHome::randomizeSpline(vector<double>& vv, const int n,
                             const vector<double>& v) {
  const double width = fabs(randNormal());
  const double height = randNormal() * v[n];  // lower limit for step

  if (n > endps.back()) {  // change gradient
    // vv[n] += GAUSS(double(n) / width) * height;
  } else {
    int ff = 0;
    double w2 = 1.0 + 3 * width;                // original verison 4.0 * width
    while (ff < nfuncs) {                       //
      if (n >= startps[ff] && n < endps[ff]) {  // update function values
        for (int i = 0; i <= w2; i++) {
          int j = n + i;
          if (j < endps[ff]) vv[j] += GAUSS(double(i) / width) * height;
          j = n - i;
          if (j >= startps[ff]) vv[j] += GAUSS(double(i) / width) * height;
        }
        break;
      }
      ff++;
    }  // while
  }    //  else
}

void pfHome::randomize(vector<double>& vv, const int n,
                       const vector<double>& v) {
  const double width = fabs(randNormal());
  const double height = randNormal() * v[n];  // lower limit for step
  double w2 = 4.0 * width;                    // original verison 4.0 * width
  for (int i = 0; i <= w2; i++) {
    int j = n + i;
    if (j < vv.size()) vv[j] += GAUSS(double(i) / width) * height;
    j = n - i;
    if (j >= 0) vv[j] += GAUSS(double(i) / width) * height;
  }
}

void pfHome::simAnneal(pfForce& fcdrv) {
  int loopcnt = 0, loopagain = 1;
  double T = dparams["temp"];
  double err = (fcdrv.*calobj[sparams["ptype"]])(ini, 1);
  double correntobj = err, overallobj = err;
  vector<double> v(nvars, dparams["istep"]);
  vector<int> naccs(nvars, 0);
  vector<double> tmpvv = encodestdv(ini);
  vector<double> optvv = encodestdv(ini);
  vector<double> errold(NEPS, err);

  while (loopcnt < iparams["kmax"] && loopagain) {
    for (int m = 0; m < NTEMP; m++) {
      for (int j = 0; j < NSTEP; j++) {
        for (int h = 0; h < nvars; h++) {
          tmpvv = ini;
          randomize(tmpvv, h, v);
          correntobj = (fcdrv.*calobj[sparams["ptype"]])(decodestdv(tmpvv), 1);

          if (correntobj <= err) {
            ini = tmpvv;
            err = correntobj;
            naccs[h]++;
            if (correntobj < overallobj) {
              optvv = tmpvv;
              overallobj = correntobj;
              (this->*write[sparams["ptype"]])();
            }
          } else if (randUniform() < (exp((err - correntobj) / T))) {
            ini = tmpvv;
            err = correntobj;
            naccs[h]++;
          }
        }  // h
      }    // steps

      for (int n = 0; n < nvars; n++) { /* step adjustment */
        if (naccs[n] > 0.6 * NSTEP)
          v[n] *= (1 + STEPVAR * ((double)naccs[n] / NSTEP - 0.6) / 0.4);
        else if (naccs[n] < 0.4 * NSTEP)
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccs[n] / NSTEP) / 0.4);
        naccs[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t\t%f\n", loopcnt, T, m + 1, error["frc"],
             overallobj);
      fflush(stdout);

      // if (iparams["resfreq"] > 0 && (m + 1) % iparams["resfreq"] == 0) {
      //   updaterho(ini);
      //   writeLMPS(optvv);
      // } else if ((m + 1) % 20 == 0)
      //   updaterho(ini);
    }  // temp

    T *= TEMPVAR;
    loopcnt++;

    for (int i = 0; i < NEPS - 1; i++) errold[i] = errold[i + 1];

    errold[NEPS - 1] = err;
    loopagain = 0;

    for (int n = 0; n < NEPS - 1; n++) {
      if (fabs(err - errold[n]) > (EPS * err * 0.01)) {
        loopagain = 1;
        break;
      }
    }

    if (!loopagain && ((err - overallobj) > (EPS * err * 0.01))) {
      ini = optvv;
      err = overallobj;
      loopagain = 1;
    }
  }  // while

  ini = optvv;
  // for growing variables
  // recorderr.push_back(overallobj);
  // recordStage(scnt++);

  v.clear();
  naccs.clear();
  tmpvv.clear();
  optvv.clear();
  errold.clear();
}

void pfHome::simAnnealSpline() {
  pfForce fcdrv(*this);
  int loopcnt = 0;
  int loopagain = 1;
  double T = dparams["temp"];
  double err = (fcdrv.*calobj[sparams["ptype"]])(decodestdv(ini), 1);
  double correntobj = err, overallobj = err;

  vector<double> v(nvars, dparams["istep"]);
  vector<int> naccs(nvars, 0);
  vector<double> tmpvv = ini;
  vector<double> optvv = ini;
  vector<double> errold(NEPS, err);

  while (loopcnt < iparams["kmax"] && loopagain) {
    for (int m = 0; m < NTEMP; m++) {
      for (int j = 0; j < NSTEP; j++) {
        for (int h = 0; h < nvars; h++) {
          tmpvv = ini;
          randomize(tmpvv, h, v);
          correntobj = (fcdrv.*calobj[sparams["ptype"]])(decodestdv(tmpvv), 1);

          if (correntobj <= err) {
            ini = tmpvv;
            err = correntobj;
            naccs[h]++;
            if (correntobj < overallobj) {
              optvv = tmpvv;
              overallobj = correntobj;
              writePot(optvv);
            }
          } else if (randUniform() < (exp((err - correntobj) / T))) {
            ini = tmpvv;
            err = correntobj;
            naccs[h]++;
          }
        }  // h
      }    // steps

      for (int n = 0; n < nvars; n++) { /* step adjustment */
        if (naccs[n] > 0.6 * NSTEP)
          v[n] *= (1 + STEPVAR * ((double)naccs[n] / NSTEP - 0.6) / 0.4);
        else if (naccs[n] < 0.4 * NSTEP)
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccs[n] / NSTEP) / 0.4);
        naccs[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\t%f\t%f\n", loopcnt, T, m + 1, error["frc"],
             error["lat"], error["ela"], overallobj);
      // error["pv"], overallobj);
      printf("%f %f %f %f\n", exprs["lat"], exprs["c11"], exprs["c12"],
             exprs["c44"]);
      // -mpcf["pv"][5].strs[0],
      // -mpcf["pv"][10].strs[0]);  // 60 GPa
      fflush(stdout);

      if (iparams["resfreq"] > 0 && (m + 1) % iparams["resfreq"] == 0) {
        updaterhoMEAM(ini, fcdrv);
        if (rescaleEMF(ini, fcdrv) == 1) {
          printf("before rescale = %f\n", err);
          err = (fcdrv.*calobj[sparams["ptype"]])(decodestdv(ini), 1);
          printf("min = %f ; max = %f ; ave = %f \n", ominrho, omaxrho,
                 oaverho);
          printf("after rescale = %f\n", err);
        }
      } else if ((m + 1) % 20 == 0)
        updaterhoMEAM(ini, fcdrv);

      // if (iparams["lmpfreq"] > 0 && (m + 1) % iparams["lmpfreq"] == 0){
      //   writeLMPS(optvv);
      //   lmpdrv->calPhy();
      // }
    }  // temp

    T *= TEMPVAR;
    loopcnt++;

    for (int i = 0; i < NEPS - 1; i++) errold[i] = errold[i + 1];

    errold[NEPS - 1] = err;
    loopagain = 0;

    for (int n = 0; n < NEPS - 1; n++) {
      if (fabs(err - errold[n]) > (EPS * err * 0.01)) {
        loopagain = 1;
        break;
      }
    }

    if (!loopagain && ((err - overallobj) > (EPS * err * 0.01))) {
      ini = optvv;
      err = overallobj;
      loopagain = 1;
    }
  }  // while

  ini = optvv;
  // recorderr.push_back(overallobj);
  // recordStage(scnt++);

  v.clear();
  naccs.clear();
  tmpvv.clear();
  optvv.clear();
  errold.clear();
}
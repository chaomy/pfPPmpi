/*
 * @Author: chaomy
 * @Date:   2017-10-30 18:46:14
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-15 22:29:08
 */

#include "pfForce.h"
#include "pfHome.h"

void pfHome::recordStage(int cnt) {
  string copy = sparams["tmpfile"];
  sparams["tmpfile"] = sparams["tmpdir"] + "/" + sparams["tmpfile"];
  sparams["tmpfile"] += ("." + to_string(cnt));
  writePot();
  sparams["tmpfile"] = copy;
  copy = sparams["lmpfile"];
  sparams["lmpfile"] = sparams["lmpdir"] + "/" + sparams["lmpfile"];
  sparams["lmpfile"] += ("." + to_string(cnt));
  writeLMPS(ini);
  sparams["lmpfile"] = copy;
}

void pfHome::writePot(const vector<double>& vv) {
  FILE* fid = fopen(sparams["tmpfile"].c_str(), "w");
  if (!fid) cerr << "error opening " + sparams["tmpfile"] << endl;
  fprintf(fid, "#F 4 %d\n", nfuncs);
  fprintf(fid, "#T %s \n", sparams["ptype"].c_str());
  fprintf(fid, "#C %s \n", sparams["elem"].c_str());
  fprintf(fid, "## %s-%s %s %s\n", sparams["elem"].c_str(),
          sparams["elem"].c_str(), sparams["elem"].c_str(),
          sparams["elem"].c_str());

  if (sparams["ptype"] == "EAM")
    fprintf(fid, "#G 3 3 3\n");
  else if (sparams["ptype"] == "ADP")
    fprintf(fid, "#G 3 3 3 3 3\n");
  else if (sparams["ptype"] == "MEAM")
    fprintf(fid, "#G 3 3 3 3 3\n");
  fprintf(fid, "#E\n");

  for (int i = 0; i < nfuncs; i++) fprintf(fid, "%d\n", funcs[i].npts);
  fprintf(fid, "\n");

  for (int i = 0; i < nfuncs; i++) {
    int stpnt = startps[i];
    fprintf(fid, "%.16e %.16e\n", funcs[i].g1.front(), funcs[i].g1.back());
    for (int j = 0; j < funcs[i].npts; j++)
      fprintf(fid, "%.16e %.16e\n", funcs[i].xx[j], vv[stpnt + j]);
    fprintf(fid, "\n");
  }
  fclose(fid);
}

void pfHome::writePot(const string& fname) {
  FILE* fid = fopen(fname.c_str(), "w");
  if (!fid) cerr << "error opening " + sparams["tmpfile"] << endl;

  fprintf(fid, "#F 4 %d\n", nfuncs);
  fprintf(fid, "#T %s \n", sparams["ptype"].c_str());
  fprintf(fid, "#C %s \n", sparams["elem"].c_str());
  fprintf(fid, "## %s-%s %s %s\n", sparams["elem"].c_str(),
          sparams["elem"].c_str(), sparams["elem"].c_str(),
          sparams["elem"].c_str());

  if (sparams["ptype"] == "EAM")
    fprintf(fid, "#G 3 3 3\n");
  else if (sparams["ptype"] == "ADP")
    fprintf(fid, "#G 3 3 3 3 3\n");
  else if (sparams["ptype"] == "MEAM")
    fprintf(fid, "#G 3 3 3 3 3\n");
  fprintf(fid, "#E\n");
  for (int i = 0; i < nfuncs; i++) fprintf(fid, "%d\n", funcs[i].npts);
  fprintf(fid, "\n");

  for (int i = 0; i < nfuncs; i++) {
    fprintf(fid, "%.16e %.16e\n", funcs[i].s.m_c0, funcs[i].s.m_c.back());
    for (int j = 0; j < funcs[i].npts; j++)
      fprintf(fid, "%.16e %.16e\n", funcs[i].xx[j], funcs[i].yy[j]);
    fprintf(fid, "\n");
  }
  fclose(fid);
}

void pfHome::writePot() {
  FILE* fid = fopen(sparams["tmpfile"].c_str(), "w");
  if (!fid) cerr << "error opening " + sparams["tmpfile"] << endl;

  fprintf(fid, "#F 4 %d\n", nfuncs);
  fprintf(fid, "#T %s \n", sparams["ptype"].c_str());
  fprintf(fid, "#C %s \n", sparams["elem"].c_str());
  fprintf(fid, "## %s-%s %s %s\n", sparams["elem"].c_str(),
          sparams["elem"].c_str(), sparams["elem"].c_str(),
          sparams["elem"].c_str());

  if (sparams["ptype"] == "EAM")
    fprintf(fid, "#G 3 3 3\n");
  else if (sparams["ptype"] == "ADP")
    fprintf(fid, "#G 3 3 3 3 3\n");
  else if (sparams["ptype"] == "MEAM")
    fprintf(fid, "#G 3 3 3 3 3\n");
  fprintf(fid, "#E\n");
  for (int i = 0; i < nfuncs; i++) fprintf(fid, "%d\n", funcs[i].npts);
  fprintf(fid, "\n");

  for (int i = 0; i < nfuncs; i++) {
    fprintf(fid, "%.16e %.16e\n", funcs[i].s.m_c0, funcs[i].s.m_c.back());
    for (int j = 0; j < funcs[i].npts; j++)
      fprintf(fid, "%.16e %.16e\n", funcs[i].xx[j], funcs[i].yy[j]);
    fprintf(fid, "\n");
  }
  fclose(fid);
}

void pfHome::writeMEAMS() {
  FILE* fid = fopen(sparams["lmpfile"].c_str(), "w");
  if (!fid) cerr << "error opening " + sparams["lmpfile"] << endl;
  fprintf(fid, "# %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", funcs[0].bl,
          funcs[0].br, funcs[1].bl, funcs[1].br, funcs[2].bl, funcs[2].br,
          funcs[3].bl, funcs[3].br, funcs[4].bl, funcs[4].br, optidx[0],
          optidx[1], optidx[2], optidx[3], optidx[4]);
  for (auto ee : smthidx) fprintf(fid, " %d", ee);
  fprintf(fid, "\n");
  fprintf(fid, "meam/spline 1 Nb\n");
  for (int i = 0; i < nfuncs; i++) {
    for (auto idx : funcs[i].rlxid) fprintf(fid, "%d ", idx);
    fprintf(fid, "\n");
    fprintf(fid, "%d\n", funcs[i].npts);
    fprintf(fid, "%.16e %.16e\n",
            funcs[i].bl == 2 ? funcs[i].s.m_c0 : funcs[i].s.m_b0,
            funcs[i].br == 2 ? funcs[i].s.m_c.back() : funcs[i].s.m_b.back());
    for (int j = 0; j < funcs[i].npts; j++)
      fprintf(fid, "%.16e %.16e %.16e\n", funcs[i].xx[j], funcs[i].yy[j],
              2. * funcs[i].s.m_b[j]);
  }
  fclose(fid);
}

void pfHome::writeLMPS(const vector<double>& vv) {
  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    Func& ff = funcs[i];
    if (i == PHI || i == RHO || i == MEAMF) {
      for (int j = 0; j < ff.npts - 1; j++) ff.yy[j] = vv[cnt++];
    } else
      for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];
    ff.s.set_points(ff.xx, ff.yy);
  }
  (!sparams["ptype"].compare("MEAM")) ? writeMEAMS() : writeLMPS();
}

void pfHome::writeLMPS() {
  FILE* fid = fopen(sparams["lmpfile"].c_str(), "w");
  pfForce fcdrv(*this);
  if (!fid) cerr << "error opening " + sparams["lmpfile"] << endl;

  // header
  fprintf(fid, "%s %s \n", sparams["elem"].c_str(), sparams["ptype"].c_str());
  fprintf(fid, "Chaoming Yang \n");
  fprintf(fid, "(%s %s) \n", __DATE__, __TIME__);
  fprintf(fid, "1 %s\n", sparams["elem"].c_str());  // n types
  /* line 5: Nrho, drho, Nr, dr, cutoff */
  double dr = funcs[PHI].xx.back() / (LMPPNTS - 1);
  double drho = funcs[EMF].xx.back() / (LMPPNTS - 1);

  fprintf(fid, "%d %lf %d %lf %lf\n", LMPPNTS, drho, LMPPNTS, dr,
          funcs[PHI].xx.back());
  fprintf(fid, "%d %lf 0\n", iparams["atomicNum"], dparams["mass"]);
  double val = 0.0;
  double r = 0.0;
  for (int i = 0; i < LMPPNTS; i++, r += drho) /* embedding function */
    fprintf(fid, "%.16e\n", funcs[EMF].s(r));

  r = 0.0;
  for (int i = 0; i < LMPPNTS; i++, r += dr) /* transfer function rho(r) */
    fprintf(fid, "%.16e\n", funcs[RHO].s(r));

  r = 0.0;
  for (int i = 0; i < LMPPNTS; i++, r += dr) /* pair function phi(r) */
    fprintf(fid, "%.16e\n", r * funcs[PHI].s(r));

  if (sparams["ptype"] == "ADP") {
    r = 0.0;
    for (int i = 0; i < LMPPNTS; i++, r += dr) { /* dipole distortion u(r) */
      fcdrv.splint(funcs[ADPU], r, val);
      fprintf(fid, "%.16e\n", val);
    }

    r = 0.0;
    for (int i = 0; i < LMPPNTS; i++, r += dr) { /* quadrupole distortion */
      fcdrv.splint(funcs[ADPW], r, val);
      fprintf(fid, "%.16e\n", val);
    }
  }  // if ADP
  fclose(fid);
}
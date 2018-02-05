/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-04 22:04:10
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

void pfHome::readMEAMC() {
  ifstream fid;
  pfUtil pfu;
  fid.open("meam.lib");
  if (!fid.is_open()) cerr << "error opening " + sparams["potfile"] << endl;
  string buff;
  vector<string> segs(1, " ");
  vector<string> tg;

  for (int i = 0; i < 3; i++) {
    getline(fid, buff);
    segs.clear();
    pfu.split(buff, " ", segs);
    for (int j = 1; j < segs.size(); j++) tg.push_back(segs[j]);
  }

  ini.clear();
  elems.clear();
  cnn1.clear();
  for (auto& ee : meamparms) ee.second.clear();
  while (getline(fid, buff)) {
    int cn = 0;
    segs.clear();
    pfu.split(buff, " ", segs);
    cout << "buff is " << buff << " " << segs.size() << endl;
    elems.push_back(segs[cn++]);
    if (!segs[cn].compare("fcc"))
      lattp.push_back(FCC);
    else if (!segs[cn].compare("bcc"))
      lattp.push_back(BCC);
    else if (!segs[cn].compare("hcp"))
      lattp.push_back(HCP);
    else if (!segs[cn].compare("dam"))
      lattp.push_back(DAM);
    else if (!segs[cn].compare("dia"))
      lattp.push_back(DIA);
    else if (!segs[cn].compare("b1"))
      lattp.push_back(B1);
    else if (!segs[cn].compare("c11"))
      lattp.push_back(C11);
    else if (!segs[cn].compare("L12"))
      lattp.push_back(L12);
    else if (!segs[cn].compare("B2"))
      lattp.push_back(B2);
    cn++;

    cnn1.push_back(stoi(segs[cn++]));
    ielement.push_back(stoi(segs[cn++]));
    atwt.push_back(stof(segs[cn++]));

    segs.clear();
    getline(fid, buff);
    pfu.split(buff, " ", segs);
    // alpha b0 b1 b2 b3 alat esub asub
    for (int i : {0, 1, 2, 3, 4, 6, 7}) ini.push_back(stof(segs[i]));

    alat.push_back(stof(segs[5]));
    segs.clear();
    getline(fid, buff);
    pfu.split(buff, " ", segs);
    //  t0, t1, t2, t3
    for (int i : {0, 1, 2, 3}) ini.push_back(stof(segs[i]));
    ini.push_back(rc_meam);  // add rc_meam to variables
    rozero.push_back(stof(segs[4]));
    ibar.push_back(stoi(segs[5]));
  }
  fid.close();
}

void pfHome::readPot() {  // read dummy.pot
  funcs.clear();
  ini.clear();

  ifstream fid;
  pfUtil pfu;

  char tmp[MAXLEN];
  fid.open(sparams["potfile"].c_str());
  if (!fid.is_open()) cerr << "error opening " + sparams["potfile"] << endl;

  string buff;
  vector<string> segs(1, " ");
  vector<int> bnds;
  vector<int> npts;

  while (getline(fid, buff)) {
    segs.clear();
    pfu.split(buff, " ", segs);
    if (!segs[0].compare("#F"))
      sscanf(buff.c_str(), "%s %s %d", tmp, tmp, &nfuncs);
    else if (!segs[0].compare("#G"))
      for (unsigned int i = 1; i < segs.size(); i++)
        bnds.push_back(stoi(segs[i]));
    else if (!segs[0].compare("#E"))
      break;
  }

  for (int i = 0; i < nfuncs; i++) {
    getline(fid, buff);
    npts.push_back(stoi(buff));
  }

  double v2[2];
  double b2[2];

  for (int i = 0; i < nfuncs; i++) {
    Func tmp;
    tmp.bnd = bnds[i];
    tmp.npts = npts[i];

    getline(fid, buff);
    getline(fid, buff);
    sscanf(buff.c_str(), "%lf %lf", &b2[0], &b2[1]);

    for (int j = 0; j < npts[i]; j++) {
      getline(fid, buff);
      sscanf(buff.c_str(), "%lf %lf", &v2[0], &v2[1]);
      tmp.xx.push_back(v2[0]);
      tmp.yy.push_back(v2[1]);
      tmp.g1.push_back(0.0);
      tmp.g2.push_back(0.0);
    }
    for (int it = 0; it < 2; it++) {
      if (b2[it] < lol.back()) b2[it] = 0.5 * lol.back();
      if (b2[it] > hil.back()) b2[it] = 0.5 * hil.back();
    }
    tmp.g1.front() = b2[0];
    tmp.g1.back() = b2[1];
    funcs.push_back(tmp);
  }
  fid.close();

  // func -> ini
  nvars = 0;
  for (Func& ff : funcs) ff.step = ff.xx[1] - ff.xx[0];
  for (int i = 0; i < nfuncs; i++) {
    Func& ff = funcs[i];
    startps.push_back(nvars);
    int nt = (i == PHI || i == RHO || i == MEAMF) ? ff.npts - 1 : ff.npts;
    for (int j = 0; j < nt; j++) {
      ini.push_back(ff.yy[j]);
      lob.push_back(lol[i]);
      hib.push_back(hil[i]);
      deb.push_back(hil[i] - lol[i]);
      nvars++;
    }
    endps.push_back(nvars);
  }  // i

  if (!sparams["ptype"].compare("MEAM"))  // boundary
    funcs[MEAMF].yy.back() = 0.0;
  funcs[PHI].yy.back() = 0.0;
  funcs[RHO].yy.back() = 0.0;
}
/*
 * @Author: chaomy
 * @Date:   2017-12-17 14:00:51
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-19 10:21:24
 */

#include "pfHome.h"

void pfHome::readMEAMS() {
  funcs.clear();
  ifstream fid;
  pfUtil pfu;
  fid.open(sparams["potfile"].c_str());
  if (!fid.is_open()) cerr << "error open " << sparams["potfile"] << endl;

  string buff;
  vector<string> segs;
  vector<int> recordbd;
  getline(fid, buff);  // read head line
  pfu.split(buff, " ", segs);

  for (int i = 1; i < 11; i++) recordbd.push_back(stoi(segs[i]));
  for (int i : {0, 1, 2, 3, 4}) optidx.push_back(stoi(segs[11 + i]));
  for (int i = 16; i < segs.size(); i++) smthidx.push_back(stoi(segs[i]));
  for (int i : {0, 1, 2, 3, 4})
    if (optidx[i] == 1) cout << "to be optimized : " << i << endl;
  for (int i : smthidx) cout << "smth idx " << i << endl;

  getline(fid, buff);  // read head line 2
  int cnt = nfuncs = 5;
  vector<int>::iterator it = recordbd.begin();
  while (--cnt >= 0) {
    segs.clear();
    Func tm;
    // find id of nodes to be relaxed
    getline(fid, buff);
    pfu.split(buff, " ", segs);
    for (auto ee : segs) tm.rlxid.push_back(stoi(ee));

    // find how many nodes in total
    getline(fid, buff);
    tm.npts = stoi(buff);
    tm.g1 = vector<double>(tm.npts, 0);
    tm.g2 = vector<double>(tm.npts, 0);
    tm.xx = vector<double>(tm.npts, 0);
    tm.yy = vector<double>(tm.npts, 0);

    // read first derivatives
    getline(fid, buff);
    sscanf(buff.c_str(), "%lf %lf", &tm.g1.front(), &tm.g1.back());

    // read xx, yy, and second derivatives
    for (int j = 0; j < tm.npts; j++) {
      getline(fid, buff);
      sscanf(buff.c_str(), "%lf %lf %lf", &tm.xx[j], &tm.yy[j], &tm.g2[j]);
    }
    tm.bl = *it++;
    tm.br = *it++;
    funcs.push_back(tm);
  }
  fid.close();
  // ricut = funcs[PHI].xx.front();
  // rocut = funcs[PHI].xx.back();
  // rhcut = funcs[RHO].xx.back();
  // for (auto ff : funcs) {
  //   cout << "npts = " << ff.npts << endl;
  //   for (int i = 0; i < ff.npts; i++)
  //     cout << ff.xx[i] << " " << ff.yy[i] << " " << ff.g2[i] << endl;
  // }
  setSplineVariables();
}
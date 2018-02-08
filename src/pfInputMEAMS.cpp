/*
 * @Author: chaomy
 * @Date:   2017-12-17 14:00:51
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-07 21:59:10
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

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

  for (int i = 1; i < segs.size(); i++) recordbd.push_back(stoi(segs[i]));

  int cnt = nfuncs = 5;
  vector<int>::iterator it = recordbd.begin();
  while (--cnt >= 0) {
    Func tm;
    getline(fid, buff);
    tm.npts = stoi(buff);
    tm.g1 = vector<double>(tm.npts, 0);
    tm.g2 = vector<double>(tm.npts, 0);
    tm.xx = vector<double>(tm.npts, 0);
    tm.yy = vector<double>(tm.npts, 0);
    getline(fid, buff);
    sscanf(buff.c_str(), "%lf %lf", &tm.g1.front(), &tm.g1.back());
    getline(fid, buff);
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

  setSplineBoundary();
}
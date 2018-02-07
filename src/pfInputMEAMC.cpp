/*
 * @Author: chaomy
 * @Date:   2018-02-06 19:10:18
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-06 19:11:22
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

void pfHome::readMEAMCcnt() {
  ifstream fid;
  pfUtil pfu;
  string buff;
  vector<string> segs;
  fid.open(sparams["meamcnt"].c_str());
  if (!fid.is_open()) cerr << "error opening " + sparams["meamcnt"] << endl;
  getline(fid, buff);
  pfu.split(buff, " ", segs);
  if (segs.size() != ini.size()) cerr << "varialbes do not match !" << endl;
  for (int i = 0; i < segs.size(); i++)
    cout << "i = " << i << " " << (ini[i] = stof(segs[i])) << endl;
  fid.close();
}

void pfHome::readMEAMC() {
  ifstream fid;
  pfUtil pfu;
  fid.open(sparams["potfile"].c_str());
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
  int in = 0;
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
    // alpha b0 b1 b2 b3 esub asub
    for (int i : {0, 1, 2, 3, 4, 6, 7}) ini.push_back(stof(segs[i]));
    alat.push_back(stof(segs[5]));
    segs.clear();
    getline(fid, buff);
    pfu.split(buff, " ", segs);
    //   t1, t2, t3  (t0 = 1)
    for (int i : {1, 2, 3}) ini.push_back(stof(segs[i]));
    t0.push_back(stof(segs[0]));
    // rc
    ini.push_back(rc_meam);
    // Cmin_meam
    ini.push_back(Cmin_meam[in][in][in]);
    // alat
    // ini.push_back(alat.back());
    rozero.push_back(stof(segs[4]));
    ibar.push_back(stoi(segs[5]));
    in++;
  }
  fid.close();
  if (!sparams["opt"].compare("cnt")) readMEAMCcnt();
}
/*
 * @Author: chaomy
 * @Date:   2018-02-06 19:10:18
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-16 16:24:57
 */

#include "pfMEAMC.h"

void pfHome::pfForce::pfMEAMC::readMEAMCcnt() {
  ifstream fid;
  string buff;
  vector<string> segs;
  fid.open(sparams["meamcnt"].c_str());
  if (!fid.is_open()) cerr << "error opening " + sparams["meamcnt"] << endl;
  getline(fid, buff);
  split(buff, " ", segs);
  if (segs.size() != ini.size()) cerr << "varialbes do not match !" << endl;
  for (int i = 0; i < segs.size(); i++)
    cout << "i = " << i << " " << (ini[i] = stof(segs[i])) << endl;
  fid.close();
}

void pfHome::pfForce::pfMEAMC::readMEAMC() {
  ifstream fid;
  fid.open(sparams["potfile"].c_str());
  if (!fid.is_open()) cerr << "error opening " + sparams["potfile"] << endl;
  string buff;
  vector<string> segs;
  vector<string> tg;

  for (int i = 0; i < 3; i++) {
    getline(fid, buff);
    segs.clear();
    split(buff, " ", segs);
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
    split(buff, " ", segs);
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
    split(buff, " ", segs);
    // alpha b0 b1 b2 b3 esub asub
    for (int i : {0, 1, 2, 3, 4, 6, 7}) ini.push_back(stof(segs[i]));
    alat.push_back(stof(segs[5]));
    segs.clear();
    getline(fid, buff);
    split(buff, " ", segs);
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

void pfHome::pfForce::pfMEAMC::writeMEAMC() {
  FILE* fid = fopen(sparams["meamlib"].c_str(), "w");
  fprintf(fid, "# elt lat z ielement atwt \n");
  fprintf(fid, "# alpha b0 b1 b2 b3 alat esub asub \n");
  fprintf(fid, "# t0 t1 t2 t3 rozero ibar\n");
  for (int i = 0; i < 1; i++) {
    fprintf(fid, "%s %s %d  %d %.4f \n", elems[i].c_str(),
            latticemp[lattp[i]].c_str(), cnn1[i], ielement[i], atwt[i]);
    fprintf(fid, "%.12f %.3f %.3f %.3f %.3f %.12f %.3f %.3f\n",
            alpha_meam[i][i], beta0_meam[i], beta1_meam[i], beta2_meam[i],
            beta3_meam[i], alat[i], Ec_meam[i][i], A_meam[i]);
    fprintf(fid, "%.1f %.2f %.2f %.2f %2.f %d\n", t0_meam[i], t1_meam[i],
            t2_meam[i], t3_meam[i], rho0_meam[i], ibar_meam[i]);
  }
  fclose(fid);

  fid = fopen("meam.param", "w");
  fprintf(fid, "rc = %f\n", rc_meam);
  fprintf(fid, "delr = %f\n", delr_meam);
  fprintf(fid, "augt1 = %d\n", augt1);
  fprintf(fid, "erose_form = %d\n", erose_form);
  fprintf(fid, "ialloy = %d\n", ialloy);
  fprintf(fid, "zbl(1,1) = %d\n", zbl_meam[0][0]);
  fprintf(fid, "nn2(1,1) = %d\n", nn2_meam[0][0]);
  fprintf(fid, "attrac(1,1) = %f\n", attrac_meam[0][0]);
  fprintf(fid, "repuls(1,1) = %f\n", repuls_meam[0][0]);
  fprintf(fid, "Cmin(1,1,1) = %f\n", Cmin_meam[0][0][0]);
  fprintf(fid, "Cmax(1,1,1) = %f\n", Cmax_meam[0][0][0]);
  fclose(fid);
}
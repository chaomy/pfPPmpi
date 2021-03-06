/*
 * @Author: chaomy
 * @Date:   2017-11-22 13:53:23
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-03 15:36:16
 */

#include "pfEle.h"

void Melem::initDFTiten() {
  /*** those stress units are kBar; remember times 0.1 to be Gpa */
  unordered_map<string, ItenT> mp;
  ItenT tmp;
  tmp.strm[0][0] = 1.000000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 1.000000;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.000000;
  tmp.egy = -20.182855;
  tmp.stsv[0] = -0.298750;
  tmp.stsv[1] = -0.298750;
  tmp.stsv[2] = -0.298750;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.000000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 1.000000;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.000000;
  tmp.egy = -40.365731;
  tmp.stsv[0] = -0.354030;
  tmp.stsv[1] = -0.459730;
  tmp.stsv[2] = -0.459730;
  tmp.stsv[3] = -0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.0] = mp;

  tmp.strm[0][0] = 1.010000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.996340;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.996800;
  tmp.egy = -20.181205;
  tmp.stsv[0] = -15.093650;
  tmp.stsv[1] = -0.465640;
  tmp.stsv[2] = -0.913380;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.010000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.995763;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.997309;
  tmp.egy = -40.362579;
  tmp.stsv[0] = -15.056340;
  tmp.stsv[1] = -0.574650;
  tmp.stsv[2] = -0.866800;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.01] = mp;

  tmp.strm[0][0] = 1.020000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.992798;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.993600;
  tmp.egy = -20.176294;
  tmp.stsv[0] = -29.814960;
  tmp.stsv[1] = -0.312350;
  tmp.stsv[2] = -1.034030;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.020000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.991550;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.994778;
  tmp.egy = -40.352808;
  tmp.stsv[0] = -29.781510;
  tmp.stsv[1] = -0.491820;
  tmp.stsv[2] = -0.902530;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.02] = mp;

  tmp.strm[0][0] = 1.030000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.989390;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.990485;
  tmp.egy = -20.168095;
  tmp.stsv[0] = -44.347500;
  tmp.stsv[1] = -0.095750;
  tmp.stsv[2] = -0.998400;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.030000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.987346;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.992620;
  tmp.egy = -40.336519;
  tmp.stsv[0] = -44.550900;
  tmp.stsv[1] = -0.682980;
  tmp.stsv[2] = -1.053890;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.03] = mp;

  tmp.strm[0][0] = 1.040000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.986101;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.987478;
  tmp.egy = -20.156760;
  tmp.stsv[0] = -58.565270;
  tmp.stsv[1] = 0.087940;
  tmp.stsv[2] = -0.950160;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.040000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.983117;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.990719;
  tmp.egy = -40.313857;
  tmp.stsv[0] = -58.968400;
  tmp.stsv[1] = -0.989070;
  tmp.stsv[2] = -1.124990;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.04] = mp;

  tmp.strm[0][0] = 1.050000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.982758;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.984798;
  tmp.egy = -20.142422;
  tmp.stsv[0] = -72.410880;
  tmp.stsv[1] = 0.220660;
  tmp.stsv[2] = -1.181660;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.050000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.978665;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.988952;
  tmp.egy = -40.285228;
  tmp.stsv[0] = -72.510290;
  tmp.stsv[1] = -0.864060;
  tmp.stsv[2] = -0.541550;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.05] = mp;

  tmp.strm[0][0] = 1.060000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.979611;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.981997;
  tmp.egy = -20.125248;
  tmp.stsv[0] = -85.563200;
  tmp.stsv[1] = 0.385860;
  tmp.stsv[2] = -1.107700;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.060000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.974013;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.987601;
  tmp.egy = -40.251050;
  tmp.stsv[0] = -85.485020;
  tmp.stsv[1] = -0.997580;
  tmp.stsv[2] = 0.023650;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.06] = mp;

  tmp.strm[0][0] = 1.070000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.976521;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.979292;
  tmp.egy = -20.105389;
  tmp.stsv[0] = -98.073050;
  tmp.stsv[1] = 0.501860;
  tmp.stsv[2] = -1.076130;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.070000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.957968;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.997638;
  tmp.egy = -40.214543;
  tmp.stsv[0] = -96.194670;
  tmp.stsv[1] = -1.131770;
  tmp.stsv[2] = 0.780260;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.07] = mp;

  tmp.strm[0][0] = 1.080000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.973465;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.976673;
  tmp.egy = -20.083010;
  tmp.stsv[0] = -109.844260;
  tmp.stsv[1] = 0.587150;
  tmp.stsv[2] = -1.068850;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.080000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.948850;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.001164;
  tmp.egy = -40.172774;
  tmp.stsv[0] = -106.310840;
  tmp.stsv[1] = -1.138280;
  tmp.stsv[2] = 0.758730;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.08] = mp;

  tmp.strm[0][0] = 1.090000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.970417;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.974132;
  tmp.egy = -20.058299;
  tmp.stsv[0] = -120.811180;
  tmp.stsv[1] = 0.654220;
  tmp.stsv[2] = -1.080990;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.090000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.940458;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.003982;
  tmp.egy = -40.126911;
  tmp.stsv[0] = -115.080600;
  tmp.stsv[1] = -1.137640;
  tmp.stsv[2] = 0.713720;
  tmp.stsv[3] = -0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.09] = mp;

  tmp.strm[0][0] = 1.100000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.967350;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.971667;
  tmp.egy = -20.031454;
  tmp.stsv[0] = -130.887170;
  tmp.stsv[1] = 0.713960;
  tmp.stsv[2] = -1.104450;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.100000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.931987;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.006811;
  tmp.egy = -40.077734;
  tmp.stsv[0] = -121.973050;
  tmp.stsv[1] = -1.179350;
  tmp.stsv[2] = 0.728530;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.10] = mp;

  tmp.strm[0][0] = 1.110000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.964228;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.969283;
  tmp.egy = -20.002690;
  tmp.stsv[0] = -140.007910;
  tmp.stsv[1] = 0.761380;
  tmp.stsv[2] = -1.135250;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.110000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.922647;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.010326;
  tmp.egy = -40.026440;
  tmp.stsv[0] = -126.110460;
  tmp.stsv[1] = -1.219000;
  tmp.stsv[2] = 0.708230;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.11] = mp;

  tmp.strm[0][0] = 1.120000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.961261;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.966732;
  tmp.egy = -19.972245;
  tmp.stsv[0] = -148.119710;
  tmp.stsv[1] = 0.720020;
  tmp.stsv[2] = -1.079850;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.120000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.911970;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.014804;
  tmp.egy = -39.974550;
  tmp.stsv[0] = -126.238030;
  tmp.stsv[1] = -1.171150;
  tmp.stsv[2] = 0.628540;
  tmp.stsv[3] = -0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.12] = mp;

  tmp.strm[0][0] = 1.130000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.957886;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.964555;
  tmp.egy = -19.940346;
  tmp.stsv[0] = -155.150560;
  tmp.stsv[1] = 0.767070;
  tmp.stsv[2] = -1.125600;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.130000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.899841;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.020099;
  tmp.egy = -39.924094;
  tmp.stsv[0] = -120.991300;
  tmp.stsv[1] = -1.166000;
  tmp.stsv[2] = 0.521480;
  tmp.stsv[3] = -0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.13] = mp;

  tmp.strm[0][0] = 1.140000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.954265;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.962553;
  tmp.egy = -19.907249;
  tmp.stsv[0] = -161.015360;
  tmp.stsv[1] = 0.803290;
  tmp.stsv[2] = -1.159300;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.140000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.886034;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.026111;
  tmp.egy = -39.877630;
  tmp.stsv[0] = -108.513690;
  tmp.stsv[1] = -1.126190;
  tmp.stsv[2] = 0.396180;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.14] = mp;

  tmp.strm[0][0] = 1.150000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.950536;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.960551;
  tmp.egy = -19.873273;
  tmp.stsv[0] = -165.642890;
  tmp.stsv[1] = 0.772490;
  tmp.stsv[2] = -1.108150;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.150000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.870554;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.032472;
  tmp.egy = -39.838482;
  tmp.stsv[0] = -86.975890;
  tmp.stsv[1] = -1.027450;
  tmp.stsv[2] = 0.280810;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.15] = mp;

  tmp.strm[0][0] = 1.160000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.945718;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.959514;
  tmp.egy = -19.838676;
  tmp.stsv[0] = -168.835780;
  tmp.stsv[1] = 0.763890;
  tmp.stsv[2] = -1.123350;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.160000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.853877;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.038563;
  tmp.egy = -39.810019;
  tmp.stsv[0] = -56.134640;
  tmp.stsv[1] = -1.112170;
  tmp.stsv[2] = 0.271710;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.16] = mp;

  tmp.strm[0][0] = 1.170000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.939152;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.960066;
  tmp.egy = -19.803814;
  tmp.stsv[0] = -170.237390;
  tmp.stsv[1] = 0.801000;
  tmp.stsv[2] = -1.150490;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.170000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.836108;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.044303;
  tmp.egy = -39.795790;
  tmp.stsv[0] = -16.410090;
  tmp.stsv[1] = -1.086420;
  tmp.stsv[2] = 0.180990;
  tmp.stsv[3] = -0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.17] = mp;

  tmp.strm[0][0] = 1.180000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.932182;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.960852;
  tmp.egy = -19.769138;
  tmp.stsv[0] = -169.774080;
  tmp.stsv[1] = 0.445170;
  tmp.stsv[2] = -1.036110;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.180000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.818491;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.049219;
  tmp.egy = -39.798286;
  tmp.stsv[0] = 26.802970;
  tmp.stsv[1] = -1.067130;
  tmp.stsv[2] = 0.101210;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.18] = mp;

  tmp.strm[0][0] = 1.190000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.925599;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.961105;
  tmp.egy = -19.735115;
  tmp.stsv[0] = -167.491670;
  tmp.stsv[1] = -0.625080;
  tmp.stsv[2] = -0.816660;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.190000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.802514;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.052961;
  tmp.egy = -39.816961;
  tmp.stsv[0] = 65.491020;
  tmp.stsv[1] = -1.103810;
  tmp.stsv[2] = 0.080860;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.19] = mp;

  tmp.strm[0][0] = 1.200000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.917681;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.961761;
  tmp.egy = -19.702359;
  tmp.stsv[0] = -161.522170;
  tmp.stsv[1] = -1.051190;
  tmp.stsv[2] = 0.534370;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.200000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.788829;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.055543;
  tmp.egy = -39.848309;
  tmp.stsv[0] = 94.981500;
  tmp.stsv[1] = -0.983210;
  tmp.stsv[2] = 0.108210;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.20] = mp;

  tmp.strm[0][0] = 1.210000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.904582;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.966898;
  tmp.egy = -19.671946;
  tmp.stsv[0] = -148.450140;
  tmp.stsv[1] = -1.031400;
  tmp.stsv[2] = 0.704500;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.210000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.777487;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.057037;
  tmp.egy = -39.888101;
  tmp.stsv[0] = 114.387280;
  tmp.stsv[1] = -0.925190;
  tmp.stsv[2] = 0.226840;
  tmp.stsv[3] = -0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.21] = mp;

  tmp.strm[0][0] = 1.220000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.887752;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 0.973562;
  tmp.egy = -19.644176;
  tmp.stsv[0] = -122.700100;
  tmp.stsv[1] = -0.953900;
  tmp.stsv[2] = 0.461720;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.220000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.768147;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.057577;
  tmp.egy = -39.932876;
  tmp.stsv[0] = 125.252780;
  tmp.stsv[1] = -0.974670;
  tmp.stsv[2] = 0.501380;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.22] = mp;

  tmp.strm[0][0] = 1.230000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.824263;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.014436;
  tmp.egy = -19.632587;
  tmp.stsv[0] = 0.993740;
  tmp.stsv[1] = -1.513790;
  tmp.stsv[2] = 0.877260;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = -0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.230000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.760160;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.057568;
  tmp.egy = -39.979952;
  tmp.stsv[0] = 129.975380;
  tmp.stsv[1] = -0.813260;
  tmp.stsv[2] = 0.631960;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.23] = mp;

  tmp.strm[0][0] = 1.240000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.816766;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.013847;
  tmp.egy = -19.635817;
  tmp.stsv[0] = 12.747850;
  tmp.stsv[1] = -1.048240;
  tmp.stsv[2] = 0.765600;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.240000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.753244;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.057187;
  tmp.egy = -40.027276;
  tmp.stsv[0] = 129.898600;
  tmp.stsv[1] = -0.529090;
  tmp.stsv[2] = 0.470440;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.24] = mp;

  tmp.strm[0][0] = 1.250000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.812251;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.010734;
  tmp.egy = -19.639651;
  tmp.stsv[0] = 18.519370;
  tmp.stsv[1] = -0.698800;
  tmp.stsv[2] = 0.615180;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["tpath"] = tmp;

  tmp.strm[0][0] = 1.250000;
  tmp.strm[0][1] = 0.000000;
  tmp.strm[0][2] = 0.000000;
  tmp.strm[1][0] = 0.000000;
  tmp.strm[1][1] = 0.747469;
  tmp.strm[1][2] = 0.000000;
  tmp.strm[2][0] = 0.000000;
  tmp.strm[2][1] = 0.000000;
  tmp.strm[2][2] = 1.056088;
  tmp.egy = -40.073375;
  tmp.stsv[0] = 125.966310;
  tmp.stsv[1] = -0.494400;
  tmp.stsv[2] = 0.692100;
  tmp.stsv[3] = 0.000000;
  tmp.stsv[4] = 0.000000;
  tmp.stsv[5] = 0.000000;
  mp["opath"] = tmp;
  itdftm[0.25] = mp;
}

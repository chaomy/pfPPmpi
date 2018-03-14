/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-03-13 16:21:38
 */

#include "pfHome.h"

void pfHome::initAngles() {
  for (Config &tmpc : configs) initAngles(tmpc);
}

void pfHome::initNeighsFull() {
  ricut = rocut;
  for (Config &tmpc : configs) {
    initBox(tmpc);
    wrapAtomPos(tmpc);
    initNeighsFull(tmpc);
  }
  if (cmm.rank() == PFROOT) cout << "rin cut = " << ricut << endl;
}

void pfHome::initBox(Config &tmpc) {
  crossProd33(tmpc.bvy, tmpc.bvz, tmpc.tvx);
  crossProd33(tmpc.bvz, tmpc.bvx, tmpc.tvy);
  crossProd33(tmpc.bvx, tmpc.bvy, tmpc.tvz);

  tmpc.vol = vecInnProd33(tmpc.bvx, tmpc.tvx);
  double inv = 1. / tmpc.vol;

  // normalize
  scaleVec(tmpc.tvx, inv);
  scaleVec(tmpc.tvy, inv);
  scaleVec(tmpc.tvz, inv);

  double iheight[3];
  iheight[0] = sqrt(square33(tmpc.tvx));
  iheight[1] = sqrt(square33(tmpc.tvy));
  iheight[2] = sqrt(square33(tmpc.tvz));

  for (int i : {X, Y, Z}) tmpc.scale[i] = (int)ceil(rocut * iheight[i]);
}

void pfHome::initAngles(Config &tmpc) {
  for (int ii = 0; ii < tmpc.natoms; ii++) {
    pfAtom &atmii = tmpc.atoms[ii];
    atmii.angMat.clear();
    for (int jj = 0; jj < atmii.nneighsFull; jj++) {
      Neigh &ngbj = atmii.neighsFull[jj];
      vector<Angle> angvec;
      for (int kk = 0; kk < jj; kk++) {
        Neigh &ngbk = atmii.neighsFull[kk];
        Angle tmpang;
        tmpang.gcos = ngbj.dist2r[X] * ngbk.dist2r[X] +
                      ngbj.dist2r[Y] * ngbk.dist2r[Y] +
                      ngbj.dist2r[Z] * ngbk.dist2r[Z];
        setAngleslotStd(tmpang, funcs[MEAMG], tmpang.gcos);
        angvec.push_back(tmpang);
      }  // kk
      atmii.angMat.push_back(angvec);
    }  // jj
  }    // ii
}

void pfHome::wrapAtomPos(Config &tmpc) {
  vector<double> boxsidelo(3), boxsidehi(3), npst(3);

  for (int k : {0, 1, 2}) {
    boxsidelo[k] = 0.0 * tmpc.bvx[k] + 0.0 * tmpc.bvy[k] + 0.0 * tmpc.bvz[k];
    boxsidehi[k] = 1.0 * tmpc.bvx[k] + 1.0 * tmpc.bvy[k] + 1.0 * tmpc.bvz[k];
  }

  for (pfAtom &atm : tmpc.atoms) {
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
        for (int iz = -1; iz <= 1; iz++) {
          if ((ix == 0) && (iy == 0) && (iz == 0)) continue;
          npst[0] = atm.pst[0] + ix * tmpc.bvx[0] + iy * tmpc.bvy[0] +
                    iz * tmpc.bvz[0];
          if (npst[0] > boxsidehi[0] || npst[0] < boxsidelo[0]) continue;

          npst[1] = atm.pst[1] + ix * tmpc.bvx[1] + iy * tmpc.bvy[1] +
                    iz * tmpc.bvz[1];
          if (npst[1] > boxsidehi[1] || npst[1] < boxsidelo[1]) continue;

          npst[2] = atm.pst[2] + ix * tmpc.bvx[2] + iy * tmpc.bvy[2] +
                    iz * tmpc.bvz[2];
          if (npst[2] > boxsidehi[2] || npst[2] < boxsidelo[2]) continue;
          // update
          for (int it : {0, 1, 2}) atm.pst[it] = npst[it];
        }
      }
    }
  }
}

void pfHome::initNeighsFull(Config &tmpc) {
  vector<double> d0(3);
  vector<double> dij(3);
  for (int ii = 0; ii < tmpc.natoms; ii++) {
    pfAtom &atmii = tmpc.atoms[ii];
    atmii.nneighsFull = 0;
    atmii.neighsFull.clear();
    int cn = 0;
    for (int jj = 0; jj < tmpc.natoms; jj++) {
      pfAtom &atmjj = tmpc.atoms[jj];

      for (int k : {0, 1, 2}) d0[k] = atmjj.pst[k] - atmii.pst[k];

      for (int ix = -tmpc.scale[0]; ix <= tmpc.scale[0]; ix++) {
        for (int iy = -tmpc.scale[1]; iy <= tmpc.scale[1]; iy++) {
          for (int iz = -tmpc.scale[2]; iz <= tmpc.scale[2]; iz++) {
            if ((ii == jj) && (ix == 0) && (iy == 0) && (iz == 0)) continue;

            for (int k : {0, 1, 2})
              dij[k] = d0[k] + ix * tmpc.bvx[k] + iy * tmpc.bvy[k] +
                       iz * tmpc.bvz[k];

            double r = sqrt(square33(dij));
            if (cmm.rank() == PFROOT && r < 2.01) cout << tmpc.cfgid << endl;

            if (r < rocut) {
              ricut = fmin(ricut, r);

              Neigh tmpn(cn++);
              double invr = 1. / r;

              tmpn.r = r;
              tmpn.r2 = r * r;
              tmpn.invr = invr;
              tmpn.aid = atmjj.id;

              tmpn.dist[X] = dij[X];
              tmpn.dist[Y] = dij[Y];
              tmpn.dist[Z] = dij[Z];

              tmpn.dist2r[X] = dij[X] * invr;
              tmpn.dist2r[Y] = dij[Y] * invr;
              tmpn.dist2r[Z] = dij[Z] * invr;

              if (sparams["ptype"] != "MEAMC") {
                setNeighslotStd(tmpn, funcs[PHI], r);
                setNeighslotStd(tmpn, funcs[RHO], r);
                setNeighslotStd(tmpn, funcs[MEAMF], r);
              }

              atmii.neighsFull.push_back(tmpn);

              if ((ix == 0) && (iy == 0) && (iz == 0)) {
                if (tmpn.aid < ii) continue;
              } else {
                if (dij[Z] < 0.0) continue;
                if (dij[Z] == 0.0) {
                  if (dij[Y] < 0.0) continue;
                  if (dij[Y] == 0.0 && dij[X] < 0.0) continue;
                }
              }
              atmii.neighidxHalf.push_back(tmpn.nid);
            }  // rcut
          }    // iz
        }      // iy
      }        // ix
    }          // jj
    atmii.nneighsFull = atmii.neighsFull.size();
  }  // ii
}

void pfHome::setNeighslotStd(Neigh &refn, Func func, double r) {
  vector<double>::const_iterator it =
      std::lower_bound(func.xx.begin(), func.xx.end(), r);
  int idx = std::max(int(it - func.xx.begin()) - 1, 0);
  refn.shifts.push_back(r - func.xx[idx]);
  if (r < func.xx.front())
    refn.slots.push_back(-1);
  else if (r > func.xx.back())
    refn.slots.push_back(-2);
  else
    refn.slots.push_back(idx);
}

void pfHome::setAngleslotStd(Angle &refang, Func func, double r) {
  vector<double>::const_iterator it =
      std::lower_bound(func.xx.begin(), func.xx.end(), r);
  int idx = std::max(int(it - func.xx.begin()) - 1, 0);
  refang.shift = (r - func.xx[idx]);
  if (r < func.xx.front())
    refang.slot = -1;
  else if (r > func.xx.back())
    refang.slot = -2;
  else
    refang.slot = idx;
}

void pfHome::setNeighslot(Neigh &refn, Func func, double r) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  int slot = 0;

  while ((hi - lo) > 1) {
    slot = (hi + lo) >> 1;
    if (func.xx[slot] > r)
      hi = slot;
    else
      lo = slot;
  }

  slot = lo;

  double step = func.xx[hi] - func.xx[lo];
  double shift = (r - func.xx[lo]) / step;

  if (slot >= (func.npts - 1)) {
    slot--;
    shift += 1;
  }

  refn.slots.push_back(slot);
  refn.steps.push_back(step);
  refn.shifts.push_back(shift);
}

void pfHome::setAngleslot(Angle &refang, Func func, double r) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  int slot = 0;

  while ((hi - lo) > 1) {
    slot = (hi + lo) >> 1;
    if (func.xx[slot] > r)
      hi = slot;
    else
      lo = slot;
  }

  slot = lo;
  double step = func.xx[hi] - func.xx[lo];
  double shift = (r - func.xx[lo]) / step;

  if (slot >= (func.npts - 1)) {
    slot--;
    shift += 1;
  }

  refang.slot = slot;
  refang.step = step;
  refang.shift = shift;
}

void pfHome::updateNeighslot(Neigh &refn, Func func, double r, int id) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  int slot = 0;

  while ((hi - lo) > 1) {
    slot = (hi + lo) >> 1;
    if (func.xx[slot] > r)
      hi = slot;
    else
      lo = slot;
  }

  slot = lo;
  double step = func.xx[hi] - func.xx[lo];
  double shift = (r - func.xx[lo]) / step;

  if (slot >= (func.npts - 1)) {
    slot--;
    shift += 1;
  }

  refn.slots[id] = slot;
  refn.steps[id] = step;
  refn.shifts[id] = shift;
}
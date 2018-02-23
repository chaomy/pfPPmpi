/*
 * @Xuthor: chaomy
 * @Date:   2018-01-10 20:08:18
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-02-22 23:52:36
 *
 * Modified from mlpack
 * Implementation of the Covariance Matrix Adaptation Evolution Strategy as
 * proposed by N. Hansen et al. in "Completely Derandomized Self-Adaptation in
 * Evolution Strategies".
 *
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

using arma::accu;
using arma::linspace;
using arma::mat;
using arma::randu;
using arma::vec;
using std::cout;
using std::endl;
using std::ofstream;
using std::to_string;

double pfHome::testFunc(arma::mat& vc) {
  return pow((vc[0] - 2.0), 2) + pow((vc[1] - 3.0), 2) + pow((vc[2] - 0.01), 2);
}

void pfHome::cntcmaes() {
  arma::mat iterate =
      5.0 + dparams["istep"] *
                (arma::mat(nvars, 1, arma::fill::randu) - 0.5);  // to [0, 10]
  (this->*calobj[sparams["ptype"]])(decodev(iterate), 1);
  if (cmm.rank() == PFROOT) {
    cmaes(iterate);
    (this->*calobj[sparams["ptype"]])(decodev(iterate), EXT);
  }
}

void pfHome::loopcmaes() {
  // start from scratch
  // arma::mat iterate(nvars, 1, arma::fill::randu);
  // iterate *= 10;
  double cr = 1e30, op = 1e30;
  arma::mat iterate =
      5.0 + dparams["istep"] *
                (arma::mat(nvars, 1, arma::fill::randu) - 0.5);  // to [0, 10]
  (this->*calobj[sparams["ptype"]])(decodev(iterate), 1);
  for (int i = 0; i < iparams["kmax"]; i++) {
    if (cmm.rank() == PFROOT) {
      phyweigh = dparams["pweight"];
      if ((cr = cmaes(iterate)) < op) {
        op = cr;
        std::rename("meam.lib.best", "meam.lib.best.overall");
        std::rename("err.txt", "err.overall");
      } else {
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << i;
        std::rename("meam.lib.best", ("meam.lib." + ss.str()).c_str());
        std::rename("err.txt", ("err." + ss.str()).c_str());
      }
      iterate = iterate + dparams["istep"] *
                              (arma::mat(nvars, 1, arma::fill::randu) - 0.5);
    }
  }
  (this->*calobj[sparams["ptype"]])(decodev(iterate), EXT);
}

double pfHome::cmaes(arma::mat& iterate) {
  int maxIt = iparams["maxstep"], lastid = 1;
  double tolerance = dparams["ftol"];
  double sigmatol = dparams["xtol"];
  ofstream of1("err.txt", std::ofstream::out);

  // Population size.
  int lambda = (4 + std::round(3 * std::log(iterate.n_elem))) * 10;
  arma::mat best(iterate);

  // Parent weights.
  const size_t mu = std::round(lambda / 2);

  cout << "nrows: " << iterate.n_rows << " lambda: " << lambda << " mu: " << mu
       << endl;

  arma::vec w = std::log(mu + 0.5) -
                arma::log(arma::linspace<arma::vec>(0, mu - 1, mu) + 1.0);
  w /= arma::sum(w);

  // Number of effective solutions.
  const double muEffective = 1 / arma::accu(arma::pow(w, 2));

  // Step size control parameters.
  arma::vec sigma(3);

  // double upperBound = 10., lowerBound = -10.;
  // sigma(0) = 0.1 * (upperBound - lowerBound);
  sigma(0) = 3. * dparams["istep"];

  const double cs = (muEffective + 2) / (iterate.n_elem + muEffective + 5);
  const double ds =
      1 + cs +
      2 * std::max(std::sqrt((muEffective - 1) / (iterate.n_elem + 1)) - 1,
                   0.0);
  const double enn =
      std::sqrt(iterate.n_elem) * (1.0 - 1.0 / (4.0 * iterate.n_elem) +
                                   1.0 / (21 * std::pow(iterate.n_elem, 2)));

  // Covariance update parameters Cumulation for distribution.
  const double cc = (4 + muEffective / iterate.n_elem) /
                    (4 + iterate.n_elem + 2 * muEffective / iterate.n_elem);
  const double h = (1.4 + 2.0 / (iterate.n_elem + 1.0)) * enn;

  const double c1 = 2 / (std::pow(iterate.n_elem + 1.3, 2) + muEffective);
  const double alphaMu = 2;
  const double cmu = std::min(
      1 - c1,
      alphaMu * (muEffective - 2 + 1 / muEffective) /
          (std::pow(iterate.n_elem + 2, 2) + alphaMu * muEffective / 2));

  arma::cube mps(iterate.n_rows, iterate.n_cols, 3);  // meam
  // mps.slice(0) =
  //     lowerBound +
  //     arma::randu(iterate.n_rows, iterate.n_cols) * (upperBound -
  //     lowerBound);
  mps.slice(0) = iterate;

  arma::mat step = arma::zeros(iterate.n_rows, iterate.n_cols);

  // Calculate the first objective function.
  double currentobj =
      (this->*calobj[sparams["ptype"]])(decodev(mps.slice(0)), 1);
  double overallobj = currentobj;
  vector<double> bestphy(4, 1e30);
  vector<string> bestfiles({"meam.lib.best", "lib.01", "lib.02", "lib.03"});
  double lastobj = 1e30;

  // Population parameters.
  arma::cube pStep(iterate.n_rows, iterate.n_cols, lambda);
  arma::cube pps(iterate.n_rows, iterate.n_cols, lambda);
  arma::vec pobj(lambda);
  arma::cube ps = arma::zeros(iterate.n_rows, iterate.n_cols, 2);
  arma::cube pc = ps;
  arma::cube C(iterate.n_elem, iterate.n_elem, 2);
  C.slice(0).eye();

  // Covariance matrix parameters.
  arma::vec eigval;
  arma::mat eigvec;
  arma::vec eigvalZero = arma::zeros(iterate.n_elem);

  // The current visitation order (sorted by population objectives).
  arma::uvec idx = arma::linspace<arma::uvec>(0, lambda - 1, lambda);

  for (size_t i = 1; i < maxIt; ++i) {
    const size_t idx0 = (i - 1) % 2;
    const size_t idx1 = i % 2;

    const arma::mat covLower = arma::chol(C.slice(idx0), "lower");

    for (size_t j = 0; j < lambda; ++j) {
      if (iterate.n_rows > iterate.n_cols) {
        pStep.slice(idx(j)) =
            covLower * arma::randn(iterate.n_rows, iterate.n_cols);
      } else {
        pStep.slice(idx(j)) =
            arma::randn(iterate.n_rows, iterate.n_cols) * covLower;
      }
      pps.slice(idx(j)) = mps.slice(idx0) + sigma(idx0) * pStep.slice(idx(j));
      pobj(idx(j)) =
          (this->*calobj[sparams["ptype"]])(decodev(pps.slice(idx(j))), 1);
    }

    // Sort population.
    idx = sort_index(pobj);

    step = w(0) * pStep.slice(idx(0));
    for (size_t j = 1; j < mu; ++j) step += w(j) * pStep.slice(idx(j));

    mps.slice(idx1) = mps.slice(idx0) + sigma(idx0) * step;

    currentobj = (this->*calobj[sparams["ptype"]])(decodev(mps.slice(idx1)), 1);

    if (currentobj < overallobj) {  // Update best parameters.
      overallobj = currentobj;
      iterate = mps.slice(idx1);

      if (i % iparams["lmpfreq"] == 0) {
        (this->*write[sparams["ptype"]])();
        lmpCheck(i, of1);
        for (int kk = 0; kk < bestphy.size(); kk++) {
          if (error["phy"] < bestphy[kk]) {
            int tmpid = kk;
            for (kk = bestphy.size() - 1; kk > tmpid; kk--) {
              bestphy[kk] = bestphy[kk - 1];
              std::rename(bestfiles[kk - 1].c_str(), bestfiles[kk].c_str());
            }
            bestphy[tmpid] = error["phy"];
            std::rename(sparams["lmpfile"].c_str(), bestfiles[tmpid].c_str());
            if (tmpid == 0) best = iterate;
            break;
          }
        }
      }
      lastid = i;
    }

    // for meams
    cout << "CMA-ES: i = " << i << ", objective " << overallobj << " "
         << error["frc"] << " " << error["engy"] << " " << error["phy"]
         << " cs " << sigma(idx1) << " " << ominrho << " " << omaxrho << " "
         << funcs[EMF].xx.front() << " " << funcs[EMF].xx.back() << " "
         << " " << configs[0].fitengy << " " << configs[0].engy << endl;

    if (iterate.n_rows > iterate.n_cols) {  // Update Step Size.
      ps.slice(idx1) =
          (1 - cs) * ps.slice(idx0) +
          std::sqrt(cs * (2 - cs) * muEffective) * covLower.t() * step;
    } else {
      ps.slice(idx1) =
          (1 - cs) * ps.slice(idx0) +
          std::sqrt(cs * (2 - cs) * muEffective) * step * covLower.t();
    }

    const double psNorm = arma::norm(ps.slice(idx1));
    sigma(idx1) =
        sigma(idx0) * std::pow(std::exp(cs / ds * psNorm / enn - 1), 0.3);

    // Update covariance matrix.
    if ((psNorm / sqrt(1 - std::pow(1 - cs, 2 * i))) < h) {
      pc.slice(idx1) = (1 - cc) * pc.slice(idx0) +
                       std::sqrt(cc * (2 - cc) * muEffective) * step;

      if (iterate.n_rows > iterate.n_cols) {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1) * pc.slice(idx1).t());
      } else {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1).t() * pc.slice(idx1));
      }
    } else {
      pc.slice(idx1) = (1 - cc) * pc.slice(idx0);

      if (iterate.n_rows > iterate.n_cols) {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1) * pc.slice(idx1).t() +
                              (cc * (2 - cc)) * C.slice(idx0));
      } else {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1).t() * pc.slice(idx1) +
                              (cc * (2 - cc)) * C.slice(idx0));
      }
    }

    if (iterate.n_rows > iterate.n_cols) {
      for (size_t j = 0; j < mu; ++j) {
        C.slice(idx1) = C.slice(idx1) + cmu * w(j) * pStep.slice(idx(j)) *
                                            pStep.slice(idx(j)).t();
      }
    } else {
      for (size_t j = 0; j < mu; ++j) {
        C.slice(idx1) = C.slice(idx1) + cmu * w(j) * pStep.slice(idx(j)).t() *
                                            pStep.slice(idx(j));
      }
    }

    arma::eig_sym(eigval, eigvec, C.slice(idx1));
    const arma::uvec negativeEigval = find(eigval < 0, 1);
    if (!negativeEigval.is_empty()) {
      if (negativeEigval(0) == 0) {
        C.slice(idx1).zeros();
      } else {
        C.slice(idx1) = eigvec.cols(0, negativeEigval(0) - 1) *
                        arma::diagmat(eigval.subvec(0, negativeEigval(0) - 1)) *
                        eigvec.cols(0, negativeEigval(0) - 1).t();
      }
    }

    if (std::isnan(overallobj) || std::isinf(overallobj)) {
      cout << "CMA-ES: converged to " << overallobj << "; "
           << "terminating with failure.  Try a smaller step size?"
           << std::endl;
      return overallobj;
    }

    if ((std::abs(lastobj - overallobj) < tolerance ||
         sigma(idx1) < sigmatol) &&
        i > 100) {
      cout << "CMA-ES: minimized within tolerance " << tolerance << "; "
           << "terminating optimization." << std::endl;
      return overallobj;
    }

    if (i % iparams["resfreq"] == 1) checkupdateBoundary(iterate);

    // if ((i % iparams["resfreq"] == 1) && checkBoundary(iterate)) {
    // cout << "update boundary" << endl;
    // updateBoundary(decodev(iterate));
    // mps.slice(idx1) =
    //     5.0 + cs * (arma::mat(nvars, 1, arma::fill::randu) - 0.5);
    // }

    lastobj = overallobj;
  }
  of1.close();

  iterate = best;  // give that has best phy constants to iterate
  return overallobj;
}

void pfHome::lmpCheck(int i, ofstream& of1) {
  error["phy"] = error["gsf"] = 0.0;
  if (!iparams["runlmp"]) {
    lmpdrv->calLatticeBCC();
    lmpdrv->calLatticeFCC();
    lmpdrv->calLatticeHCP();
    lmpdrv->calGSFUrlx();
  }
  lmpdrv->calElastic();
  lmpdrv->calSurface();

  remove("no");
  remove("log.lammps");
  remove("restart.equil");
  lmpdrv->exprs["bcc2hcp"] = lmpdrv->exprs["ehcp"] - lmpdrv->exprs["ebcc"];
  lmpdrv->exprs["bcc2fcc"] = lmpdrv->exprs["efcc"] - lmpdrv->exprs["ebcc"];
  vector<string> aa({"lat", "c11", "c12", "c44", "suf110", "suf100", "suf111",
                     "bcc2fcc", "bcc2hcp"});

  vector<double> ww({1e5, 400, 1000, 3000, 1e3, 1e3, 1e3, 5e3, 5e3});
  for (int i = 0; i < aa.size(); i++) {
    string ee(aa[i]);
    error["phy"] +=
        (lmpdrv->error[ee] =
             ww[i] * square11((lmpdrv->exprs[ee] - lmpdrv->targs[ee]) /
                              lmpdrv->targs[ee]));
  }

  for (int i : lmpdrv->gsfpnts)
    error["gsf"] +=
        50 * (lmpdrv->lgsf["111e110"][i] + lmpdrv->lgsf["111e211"][i]);
  error["phy"] += error["gsf"];

  of1 << i << "   " << std::setprecision(4) << lmpdrv->exprs["lat"] << " "
      << lmpdrv->exprs["c11"] << " " << lmpdrv->exprs["c12"] << " "
      << lmpdrv->exprs["c44"] << " " << lmpdrv->exprs["suf110"] << " "
      << lmpdrv->exprs["suf100"] << " " << lmpdrv->exprs["suf111"] << " "
      << lmpdrv->exprs["bcc2fcc"] << " " << lmpdrv->exprs["bcc2hcp"] << " "
      << lmpdrv->lgsf["111z110"][5] << " " << lmpdrv->lgsf["111z211"][5] << "  "
      << error["phy"] << " " << lmpdrv->error["lat"] << " "
      << lmpdrv->error["c11"] << " " << lmpdrv->error["c12"] << " "
      << lmpdrv->error["c44"] << " " << lmpdrv->error["suf110"] << " "
      << lmpdrv->error["suf100"] << " " << lmpdrv->error["suf111"] << " "
      << lmpdrv->error["bcc2fcc"] << " " << lmpdrv->error["bcc2hcp"] << " "
      << error["gsf"] << endl;
}

// for meamc
// if (i == 1) {
// of2 << i << " " << std::setprecision(4) << alpha_meam[0][0] << " "
//     << beta0_meam[0] << " " << beta1_meam[0] << " " << beta2_meam[0] << "
//     "
//     << beta3_meam[0] << " " << Ec_meam[0][0] << " " << A_meam[0] << " "
//     << t1_meam[0] << " " << t2_meam[0] << " " << t3_meam[0] << " "
//     << rc_meam << " " << Cmin_meam[0][0][0] << " "
//     << re_meam[0][0] * 2. / sqrt(3) << endl;
// }
#ifndef LIBKRIGING_SRC_LIB_INCLUDE_LIBKRIGING_LINEARALGEBRA_HPP
#define LIBKRIGING_SRC_LIB_INCLUDE_LIBKRIGING_LINEARALGEBRA_HPP

#include "libKriging/utils/lk_armadillo.hpp"

#include "libKriging/libKriging_exports.h"

class LinearAlgebra {
 public:
  static const arma::solve_opts::opts default_solve_opts;

  static double num_nugget;
  LIBKRIGING_EXPORT static void set_num_nugget(double nugget);
  LIBKRIGING_EXPORT static double get_num_nugget();

  static bool warn_chol;
  LIBKRIGING_EXPORT static void set_chol_warning(bool warn);

  static int max_inc_choldiag;
  LIBKRIGING_EXPORT static arma::mat safe_chol_lower(arma::mat X);
  static arma::mat safe_chol_lower(arma::mat X, int warn);

  static double min_rcond;
  LIBKRIGING_EXPORT static double rcond_chol(arma::mat chol);
  static double min_rcond_approx;
  LIBKRIGING_EXPORT static double rcond_approx_chol(arma::mat chol);

  LIBKRIGING_EXPORT static arma::mat cholCov(arma::mat* R,
                                      const arma::mat& _dX,
                                      const arma::vec& _theta,
                                      std::function<double(const arma::vec&, const arma::vec&)> Cov);
  LIBKRIGING_EXPORT static arma::mat update_cholCov(arma::mat* R,
                                        const arma::mat& _dX,
                                        const arma::vec& _theta, 
                                        std::function<double(const arma::vec&, const arma::vec&)> Cov,
                                        const arma::mat& Told);

  LIBKRIGING_EXPORT static arma::mat chol_block(const arma::mat C, const arma::mat Loo, const arma::mat Coo);
};

#endif  // LIBKRIGING_SRC_LIB_INCLUDE_LIBKRIGING_LINEARALGEBRA_HPP

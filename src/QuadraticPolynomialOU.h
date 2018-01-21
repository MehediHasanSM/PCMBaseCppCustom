/*
 *  QuadraticPolynomialOU.h
 *  PCMBaseCpp
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of PCMBaseCpp: A C++ backend for calculating the likelihood of phylogenetic comparative models.
 *
 * PCMBaseCpp is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU General Public License as
 * published by the Free Software Foundation.
 *
 * PCMBaseCpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with PCMBaseCpp.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */
#ifndef QuadraticPolynomial_OU_H_
#define QuadraticPolynomial_OU_H_

#include "QuadraticPolynomial.h"
#include <armadillo>
#include <sstream>

namespace PCMBaseCpp {


typedef splittree::OrderedTree<splittree::uint, LengthAndRegime> OUTreeType;

class OU: public QuadraticPolynomial<OUTreeType> {
public:
  typedef OUTreeType TreeType;
  typedef QuadraticPolynomial<TreeType> BaseType;
  typedef OU MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef splittree::PostOrderTraversal<MyType> AlgorithmType;



  // a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
  //  eigenvalues of the parameter matrix Alpha. This threshold-values is used as a condition to
  // take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
  //   `(Lambda_i+Lambda_j) --> 0`.
  double threshold_Lambda_ij_ = 1e-8;
  
  
  //
  // model parameters
  //
  
  // number of regimes;
  uint R; 
  
  
  // Each slice or column of the following cubes or matrices correponds to one regime
  arma::cube Alpha;
  arma::mat Theta;
  arma::cube Sigma;
  arma::cube Sigmae;

  arma::cx_cube P;
  arma::cx_cube P_1;
  arma::cx_cube P_1SigmaP_1_t;

  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda;

  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime
  arma::cx_cube Lambda_ij;

  arma::cube e_At;
  arma::mat I;

  OU(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data) {

    arma::uword k = BaseType::k;
    this->I = arma::eye(k, k);
  }

  void SetParameter(ParameterType const& par) {
    // The number of regimes R is deduced from the length of the par vector and
    // the number of traits, BaseType::k.
    using namespace arma;
    uword k = BaseType::k;

    if(par.size() % (k*k + k + k*k + k*k) != 0) {
      std::ostringstream os;
      os<<"ERR:03301:PCMBaseCpp:QuadraticPolynomialOU.h:SetParameter:: The length of the parameter vector ("<<par.size()<<
        ") should be a multiple of 3k^2+k, where k="<<k<<
          " is the number of traits.";
      throw std::logic_error(os.str());
    }

    this->R = par.size() / (k*k + k + k*k + k*k);

    this->Alpha = cube(&par[0], k, k, R);
    this->Theta = mat(&par[k*k*R], k, R);
    this->Sigma = cube(&par[k*k*R + k*R], k, k, R);
    this->Sigmae = cube(&par[k*k*R + k*R + k*k*R], k, k, R);

    this->P = cx_cube(k, k, R);
    this->P_1 = cx_cube(k, k, R);
    this->P_1SigmaP_1_t = cx_cube(k, k, R);
    this->lambda = cx_mat(k, R);
    this->Lambda_ij = cx_cube(k, k, R);

    this->e_At = cube(k, k, this->ref_tree_.num_nodes());

    for(uword r = 0; r < R; ++r) {
      using namespace std;

      cx_vec eigval;
      cx_mat eigvec;

      eig_gen(eigval, eigvec, Alpha.slice(r));

      lambda.col(r) = eigval;
      P.slice(r) = eigvec;

      P_1.slice(r) = inv(P.slice(r));

      P_1SigmaP_1_t.slice(r) = P_1.slice(r) * Sigma.slice(r) * P_1.slice(r).t();

      for(uword i = 0; i < k; ++i)
        for(uword j = i; j < k; ++j) {
          Lambda_ij.slice(r)(i,j) = Lambda_ij.slice(r)(j,i) =
            lambda.col(r)(i) + lambda.col(r)(j);
        }
    }
  }

  inline void InitNode(splittree::uint i) {
    BaseType::InitNode(i);

    using namespace arma;

    if(i < this->ref_tree_.num_nodes() - 1) {
      uword k = BaseType::k;

      splittree::uint j = this->ref_tree_.FindIdOfParent(i);

      double ti = this->ref_tree_.LengthOfBranch(i).length_;
      uword ri = this->ref_tree_.LengthOfBranch(i).regime_;

      uvec kj = BaseType::pc[j];
      uvec ki = BaseType::pc[i];

      arma::cx_mat fLambda_ij(k, k);

      for(uword w = 0; w < k; ++w)
        for(uword j = w; j < k; ++j) {
          if(abs(Lambda_ij.slice(ri)(w,j)) <= threshold_Lambda_ij_) {
            fLambda_ij(w,j) = fLambda_ij(j,w) = ti;
          } else {
            fLambda_ij(w,j) = fLambda_ij(j,w) =
              (1.0 - exp(-Lambda_ij.slice(ri)(w,j) * ti)) / Lambda_ij.slice(ri)(w,j);
          }
        }

      V.slice(i) = real(P.slice(ri) * (fLambda_ij % P_1SigmaP_1_t.slice(ri)) * P.slice(ri).t());

      if(i < this->ref_tree_.num_tips()) {
        V.slice(i) += Sigmae.slice(ri);
      }

      arma::uvec ui(1);
      ui(0) = i;

      V_1.slice(i)(ki, ki) = inv(V.slice(i)(ki,ki));
      //e_At.slice(i) = expmat(-ti*Alpha.slice(ri));
      e_At.slice(i) = real(P.slice(ri) * diagmat(exp(-ti * lambda.col(ri))) * P_1.slice(ri));

      A.slice(i)(ki,ki) = -0.5*V_1.slice(i)(ki,ki);
      b(ki,ui) = V_1.slice(i)(ki,ki) * (I.rows(ki) - e_At.slice(i).rows(ki)) * Theta.col(ri);
      C.slice(i)(kj,kj) = -0.5*e_At.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki) * e_At.slice(i)(ki,kj);
      d(kj,ui) = -e_At.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki) * (I.rows(ki)-e_At.slice(i).rows(ki)) * Theta.col(ri);
      E.slice(i)(kj,ki) = e_At.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki);
      f(i) =
        -0.5*(ki.n_elem * M_LN_2PI + log(det(V.slice(i)(ki,ki))) +
        Theta.col(ri).t() * (I.rows(ki)-e_At.slice(i).rows(ki)).t() *
        V_1.slice(i)(ki,ki) * (I.rows(ki)-e_At.slice(i).rows(ki)) * Theta.col(ri)).at(0,0);
    }
  }
};


typedef splittree::TraversalTask<OU> QuadraticPolynomialOU;
}

#endif // QuadraticPolynomial_OU_H_

/*
 *  QuadraticPolynomialJOU.h
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
#ifndef QuadraticPolynomial_JOU_H_
#define QuadraticPolynomial_JOU_H_

#include "QuadraticPolynomial.h"
#include <armadillo>
#include <sstream>
#include <iostream>

namespace PCMBaseCpp {


typedef splittree::OrderedTree<splittree::uint, LengthRegimeAndJump> JOUTreeType;

class JOU: public QuadraticPolynomial<JOUTreeType> {
public:
  typedef JOUTreeType TreeType;
  typedef QuadraticPolynomial<TreeType> BaseType;
  typedef JOU MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef splittree::PostOrderTraversal<MyType> AlgorithmType;

  // A 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
  // eigenvalues of the parameter matrix H. This threshold-values is used as a condition to
  // take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
  // `(Lambda_i+Lambda_j) --> 0`.
  //
  double threshold_Lambda_ij_ = 1e-8;
  
  //
  // model parameters
  //
  
  // number of regimes;
  uint R; 
    
  // Each slice or column of the following cubes or matrices correponds to one regime
  arma::cube H;
  arma::mat Theta;
  arma::cube Sigma;
  arma::cube Sigmae;

  // Jump mean and standard variance covariance matrix
  arma::mat mj;
  arma::cube Sigmaj;
  
  arma::cx_cube P;
  arma::cx_cube P_1;
  arma::cx_cube P_1SigmaP_1_t;

  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda;

  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime
  arma::cx_cube Lambda_ij;

  arma::cube e_Ht;
  arma::mat I;

  JOU(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data) {

    arma::uword k = BaseType::k;
    this->I = arma::eye(k, k);
  }

  void SetParameter(ParameterType const& par) {
    // The number of regimes R is deduced from the length of the par vector and
    // the number of traits, BaseType::k.
    using namespace arma;
    uword k = BaseType::k;

    if(par.size() % (k*k + k + k*k + k*k + k + k*k) != 0) {
      std::ostringstream os;
      os<<"ERR:03401:PCMBaseCpp:QuadraticPolynomialJOU.h:SetParameter:: The length of the parameter vector ("<<par.size()<<
        ") should be a multiple of 4k^2+2k, where k="<<k<<
          " is the number of traits.";
      throw std::logic_error(os.str());
    }

    
    this->R = par.size() / (k*k + k + k*k + k*k + k + k*k);

    this->H = cube(&par[0], k, k, R);
    this->Theta = mat(&par[k*k*R], k, R);
    this->Sigma = cube(&par[k*k*R + k*R], k, k, R);
    this->Sigmae = cube(&par[k*k*R + k*R + k*k*R], k, k, R);

    this->mj = mat(&par[k*k*R + k*R + k*k*R + k*k*R], k, R);
    this->Sigmaj = cube(&par[k*k*R + k*R + k*k*R + k*k*R + k*R], k, k, R);
    
    this->P = cx_cube(k, k, R);
    this->P_1 = cx_cube(k, k, R);
    this->P_1SigmaP_1_t = cx_cube(k, k, R);
    this->lambda = cx_mat(k, R);
    this->Lambda_ij = cx_cube(k, k, R);

    this->e_Ht = cube(k, k, this->ref_tree_.num_nodes());

    for(uword r = 0; r < R; ++r) {
      using namespace std;

      cx_vec eigval;
      cx_mat eigvec;

      eig_gen(eigval, eigvec, H.slice(r));

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
      u8 xi = this->ref_tree_.LengthOfBranch(i).jump_;

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

      //e_Ht.slice(i) = expmat(-ti*H.slice(ri));
      e_Ht.slice(i) = real(P.slice(ri) * diagmat(exp(-ti * lambda.col(ri))) * P_1.slice(ri));
      
      V.slice(i) = 
        real(P.slice(ri)* (fLambda_ij%P_1SigmaP_1_t.slice(ri)) * P.slice(ri).t() +
        xi * e_Ht.slice(i) * Sigmaj.slice(ri) * e_Ht.slice(i).t() );

      if(i < this->ref_tree_.num_tips()) {
        V.slice(i) += Sigmae.slice(ri);
      }

      V_1.slice(i)(ki, ki) = inv(V.slice(i)(ki,ki));

      arma::uvec ui(1);
      ui(0) = i;
      
      A.slice(i)(ki,ki) = -0.5*V_1.slice(i)(ki,ki);
      b(ki,ui) = V_1.slice(i)(ki,ki) * ((I.rows(ki) - e_Ht.slice(i).rows(ki)) * Theta.col(ri) + xi * (e_Ht.slice(i).rows(ki) * mj.col(ri)));
      C.slice(i)(kj,kj) = -0.5*e_Ht.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki) * e_Ht.slice(i)(ki,kj);
      d(kj,ui) = -e_Ht.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki) * ((I.rows(ki)-e_Ht.slice(i).rows(ki)) * Theta.col(ri) + xi * (e_Ht.slice(i).rows(ki) * mj.col(ri)));
      E.slice(i)(kj,ki) = e_Ht.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki);
      f(i) =
        -0.5*(ki.n_elem * M_LN_2PI + log(det(V.slice(i)(ki,ki))) +
        (Theta.col(ri).t() * (I.rows(ki)-e_Ht.slice(i).rows(ki)).t() + xi*mj.col(ri).t()*e_Ht.slice(i).rows(ki).t()) * V_1.slice(i)(ki,ki) * ((I.rows(ki)-e_Ht.slice(i).rows(ki)) * Theta.col(ri) + xi*e_Ht.slice(i).rows(ki)*mj.col(ri))).at(0,0);
    }
  }
};


typedef splittree::TraversalTask<JOU> QuadraticPolynomialJOU;
}

#endif // QuadraticPolynomial_JOU_H_

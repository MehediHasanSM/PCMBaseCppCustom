/*
 *  QuadraticPolynomialTwoSpeedOU.h
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
#ifndef QuadraticPolynomial_TwoSpeedOU_H_
#define QuadraticPolynomial_TwoSpeedOU_H_

#include "QuadraticPolynomial.h"
#include <armadillo>
#include <sstream>

namespace PCMBaseCpp {


typedef splittree::OrderedTree<splittree::uint, LengthAndRegime> OUTreeType;

class TwoSpeedOU: public QuadraticPolynomial<OUTreeType> {
public:
  typedef OUTreeType TreeType;
  typedef QuadraticPolynomial<TreeType> BaseType;
  typedef TwoSpeedOU MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef splittree::PostOrderTraversal<MyType> AlgorithmType;



  //
  // model parameters
  //
  
  // number of regimes;
  uint R; 

  // Each slice or column of the following cubes or matrices correponds to one regime
  arma::cube Alpha1;
  arma::cube Alpha2;
  arma::mat Theta;
  arma::cube Sigma;
  arma::cube Sigmae;

  // for Alpha1 defining the attraction to Theta
  arma::cx_cube P1;
  arma::cx_cube P1_1;
  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda1;
  
  // for Alpha2 defining the rate of decorrelation
  arma::cx_cube P2;
  arma::cx_cube P2_1;
  arma::cx_cube P2_1SigmaP2_1_t;

  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda2;

  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime
  arma::cx_cube Lambda2_ij;

  arma::cube V;
  arma::cube V_1;
  arma::cube e_A1t;
  arma::mat I;

  TwoSpeedOU(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data) {

    arma::uword k = BaseType::k;
    this->I = arma::eye(k, k);
  }

  void SetParameter(ParameterType const& par) {
    // The number of regimes R is deduced from the length of the par vector and
    // the number of traits, BaseType::k.
    using namespace arma;
    uword k = BaseType::k;

    if(par.size() % (2*k*k + k + k*k + k*k) != 0) {
      std::ostringstream os;
      os<<"The length of the parameter vector ("<<par.size()<<
        ") should be a multiple of 4k^2+k, where k="<<k<<
          " is the number of traits.";
      throw std::logic_error(os.str());
    }

    this->R = par.size() / (2*k*k + k + k*k + k*k);

    this->Alpha1 = cube(&par[0], k, k, R);
    this->Alpha2 = cube(&par[k*k*R], k, k, R);
    this->Theta = mat(&par[2*k*k*R], k, R);
    this->Sigma = cube(&par[2*k*k*R + k*R], k, k, R);
    this->Sigmae = cube(&par[2*k*k*R + k*R + k*k*R], k, k, R);

    this->P1 = cx_cube(k, k, R);
    this->P1_1 = cx_cube(k, k, R);
    this->lambda1 = cx_mat(k, R);
    
    this->P2 = cx_cube(k, k, R);
    this->P2_1 = cx_cube(k, k, R);
    this->P2_1SigmaP2_1_t = cx_cube(k, k, R);
    this->lambda2 = cx_mat(k, R);
    this->Lambda2_ij = cx_cube(k, k, R);

    this->V = cube(k, k, this->ref_tree_.num_nodes());
    this->V_1 = cube(k, k, this->ref_tree_.num_nodes());
    
    this->e_A1t = cube(k, k, this->ref_tree_.num_nodes());

    for(uword r = 0; r < R; ++r) {
      using namespace std;

      cx_vec eigval1;
      cx_mat eigvec1;
      eig_gen(eigval1, eigvec1, Alpha1.slice(r));
      lambda1.col(r) = eigval1;
      P1.slice(r) = eigvec1;
      P1_1.slice(r) = inv(P1.slice(r));
      
      cx_vec eigval2;
      cx_mat eigvec2;
      eig_gen(eigval2, eigvec2, Alpha2.slice(r));
      lambda2.col(r) = eigval2;
      P2.slice(r) = eigvec2;
      P2_1.slice(r) = inv(P2.slice(r));
      P2_1SigmaP2_1_t.slice(r) = P2_1.slice(r) * Sigma.slice(r) * P2_1.slice(r).t();

      for(uword i = 0; i < k; ++i)
        for(uword j = i; j < k; ++j) {
          Lambda2_ij.slice(r)(i,j) = Lambda2_ij.slice(r)(j,i) =
            lambda2.col(r)(i) + lambda2.col(r)(j);
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

      arma::cx_mat fLambda2_ij(k, k);

      for(uword w = 0; w < k; ++w)
        for(uword j = w; j < k; ++j) {
          if(abs(Lambda2_ij.slice(ri)(w,j)) == 0) {
            fLambda2_ij(w,j) = fLambda2_ij(j,w) = ti;
          } else {
            fLambda2_ij(w,j) = fLambda2_ij(j,w) =
              (1.0 - exp(-Lambda2_ij.slice(ri)(w,j) * ti)) / Lambda2_ij.slice(ri)(w,j);
          }
        }

      V.slice(i) = real(P2.slice(ri) * (fLambda2_ij % P2_1SigmaP2_1_t.slice(ri)) * P2.slice(ri).t());

      if(i < this->ref_tree_.num_tips()) {
        V.slice(i) += Sigmae.slice(ri);
      }

      arma::uvec ui(1);
      ui(0) = i;

      V_1.slice(i)(ki, ki) = inv(V.slice(i)(ki,ki));
      //e_At.slice(i) = expmat(-ti*Alpha1.slice(ri));
      e_A1t.slice(i) = real(P1.slice(ri) * diagmat(exp(-ti * lambda1.col(ri))) * P1_1.slice(ri));

      this->A.slice(i)(ki,ki) = -0.5*V_1.slice(i)(ki,ki);
      this->b(ki,ui) = V_1.slice(i)(ki,ki) * (I.rows(ki) - e_A1t.slice(i).rows(ki)) * Theta.col(ri);
      C.slice(i)(kj,kj) = -0.5*e_A1t.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki) * e_A1t.slice(i)(ki,kj);
      d(kj,ui) = -e_A1t.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki) * (I.rows(ki)-e_A1t.slice(i).rows(ki)) * Theta.col(ri);
      E.slice(i)(kj,ki) = e_A1t.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki);
      f(i) =
        -0.5*(ki.n_elem * M_LN_2PI + log(det(V.slice(i)(ki,ki))) +
        Theta.col(ri).t() * (I.rows(ki)-e_A1t.slice(i).rows(ki)).t() *
        V_1.slice(i)(ki,ki) * (I.rows(ki)-e_A1t.slice(i).rows(ki)) * Theta.col(ri)).at(0,0);
    }
  }
};


typedef splittree::TraversalTask<TwoSpeedOU> QuadraticPolynomialTwoSpeedOU;
}

#endif // QuadraticPolynomial_OU_H_

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

template<class TreeType, class DataType>
struct CondGaussianTwoSpeedOU: public CondGaussianOmegaPhiV {
  
  TreeType const& ref_tree_;
  
  // a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
  //  eigenvalues of the parameter matrix H. This threshold-values is used as a condition to
  // take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
  //   `(Lambda_i+Lambda_j) --> 0`.
  double threshold_Lambda_ij_ = 1e-8;
  
  // number of traits
  uint k_;
  
  // number of regimes;
  uint R_;
  
  //
  // model parameters
  //
 
  
  // Each slice or column of the following cubes or matrices correponds to one regime
  arma::mat X0;
  
  arma::cube H1;
  arma::cube H2;
  arma::mat Theta;
  arma::cube Sigma;
  arma::cube Sigmae;
  
  // for H1 defining the attraction to Theta
  arma::cx_cube P1;
  arma::cx_cube P1_1;
  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda1;
  
  // for H2 defining the rate of decorrelation
  arma::cx_cube P2;
  arma::cx_cube P2_1;
  arma::cx_cube P2_1SigmaP2_1_t;
  
  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda2;
  
  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime
  arma::cx_cube Lambda2_ij;
  
  arma::cube e_H1t;
  arma::mat I;
  
  CondGaussianTwoSpeedOU(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = R;
    
    this->threshold_Lambda_ij_ = ref_data.threshold_Lambda_ij_;
    
    this->I = arma::eye(k_, k_);
  }
  
  
  CondGaussianTwoSpeedOU(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = ref_data.R_;
    
    this->threshold_Lambda_ij_ = ref_data.threshold_Lambda_ij_;
    
    this->I = arma::eye(k_, k_);
  }
  
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;
    
    uint npar = R_*(4*k_*k_ + 2*k_);
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os<<"ERR:03401:PCMBaseCpp:QuadraticPolynomialTwoSpeedOU.h:CondTwoSpeedOU.SetParameter:: The length of the parameter vector minus offset ("<<par.size() - offset<<
        ") should be at least of R*(4k^2+2k), where k="<<k_<<" is the number of traits and "<<
          " R="<<R_<<" is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    this->X0 = mat(&par[offset], k_, R_);
    this->H1 = cube(&par[offset + k_*R_], k_, k_, R_);
    this->H2 = cube(&par[offset + (k_ + k_*k_)*R_], k_, k_, R_);
    this->Theta = mat(&par[offset + (k_ + k_*k_ + k_*k_)*R_], k_, R_);
    this->Sigma = cube(&par[offset + (k_ + k_*k_ + k_*k_ + k_)*R_], k_, k_, R_);
    this->Sigmae = cube(&par[offset + (k_ + k_*k_ + k_*k_ + k_ + k_*k_)*R_], k_, k_, R_);
    
    this->P1 = cx_cube(k_, k_, R_);
    this->P1_1 = cx_cube(k_, k_, R_);
    this->lambda1 = cx_mat(k_, R_);
    
    this->P2 = cx_cube(k_, k_, R_);
    this->P2_1 = cx_cube(k_, k_, R_);
    this->P2_1SigmaP2_1_t = cx_cube(k_, k_, R_);
    this->lambda2 = cx_mat(k_, R_);
    this->Lambda2_ij = cx_cube(k_, k_, R_);
    
    this->e_H1t = cube(k_, k_, this->ref_tree_.num_nodes());
    
    for(uword r = 0; r < R_; ++r) {
      using namespace std;
      
      cx_vec eigval1;
      cx_mat eigvec1;
      eig_gen(eigval1, eigvec1, H1.slice(r));
      lambda1.col(r) = eigval1;
      P1.slice(r) = eigvec1;
      P1_1.slice(r) = inv(P1.slice(r));
      
      cx_vec eigval2;
      cx_mat eigvec2;
      eig_gen(eigval2, eigvec2, H2.slice(r));
      lambda2.col(r) = eigval2;
      P2.slice(r) = eigvec2;
      P2_1.slice(r) = inv(P2.slice(r));
      P2_1SigmaP2_1_t.slice(r) = P2_1.slice(r) * Sigma.slice(r) * P2_1.slice(r).t();
      
      for(uword i = 0; i < k_; ++i)
        for(uword j = i; j < k_; ++j) {
          Lambda2_ij.slice(r)(i,j) = Lambda2_ij.slice(r)(j,i) =
            lambda2.col(r)(i) + lambda2.col(r)(j);
        }
    }
    return npar;
  }
  
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::mat& omega, arma::cube& Phi, arma::cube& V) {
    using namespace arma;
    
    double ti = this->ref_tree_.LengthOfBranch(i).length_;
    
    arma::cx_mat fLambda2_ij(k_, k_);
    
    for(uword w = 0; w < k_; ++w)
      for(uword j = w; j < k_; ++j) {
        if(abs(Lambda2_ij.slice(ri)(w,j)) <= threshold_Lambda_ij_) {
          fLambda2_ij(w,j) = fLambda2_ij(j,w) = ti;
        } else {
          fLambda2_ij(w,j) = fLambda2_ij(j,w) =
            (1.0 - exp(-Lambda2_ij.slice(ri)(w,j) * ti)) / Lambda2_ij.slice(ri)(w,j);
        }
      }
      
      Phi.slice(i) = real(P1.slice(ri) * diagmat(exp(-ti * lambda1.col(ri))) * P1_1.slice(ri));
    omega.col(i) = (I - Phi.slice(i)) * Theta.col(ri);
    
    V.slice(i) = real(P2.slice(ri) * (fLambda2_ij % P2_1SigmaP2_1_t.slice(ri)) * P2.slice(ri).t());
    if(i < this->ref_tree_.num_tips()) {
      V.slice(i) += Sigmae.slice(ri);
    }
  }
};

class TwoSpeedOU: public QuadraticPolynomial<OUTreeType> {
public:
  typedef OUTreeType TreeType;
  typedef QuadraticPolynomial<TreeType> BaseType;
  typedef TwoSpeedOU MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef splittree::PostOrderTraversal<MyType> AlgorithmType;

  CondGaussianTwoSpeedOU<TreeType, DataType> cond_dist_;
  
  TwoSpeedOU(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }
  
  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};


typedef splittree::TraversalTask<TwoSpeedOU> QuadraticPolynomialTwoSpeedOU;
}

#endif // QuadraticPolynomial_OU_H_

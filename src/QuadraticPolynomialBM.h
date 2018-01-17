/*
 *  QuadraticPolynomialBM.h
 *  PCMBaseCpp
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of PCMBaseCpp: A C++ backend for calculating the likelihood
 *  of phylogenetic comparative models.
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
#ifndef QuadraticPolynomial_BM_H_
#define QuadraticPolynomial_BM_H_

#include "QuadraticPolynomial.h"
#include <armadillo>
#include <sstream>

namespace PCMBaseCpp {


typedef splittree::OrderedTree<splittree::uint, LengthAndRegime> BMTreeType;

class BM: public QuadraticPolynomial<BMTreeType> {
public:
  typedef BMTreeType TreeType;
  typedef QuadraticPolynomial<TreeType> BaseType;
  typedef BM MyType;
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
  arma::cube Sigma;
  arma::cube Sigmae;

  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime

  arma::cube V;
  arma::cube V_1;
  arma::mat I;

  BM(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data) {

    arma::uword k = BaseType::k;
    this->I = arma::eye(k, k);
  }

  void SetParameter(ParameterType const& par) {
    // The number of regimes R is deduced from the length of the par vector and
    // the number of traits, BaseType::k.
    using namespace arma;
    uword k = BaseType::k;

    if(par.size() % (k*k + k*k) != 0) {
      std::ostringstream os;
      os<<"The length of the parameter vector ("<<par.size()<<
        ") should be a multiple of 3k^2+k, where k="<<k<<
          " is the number of traits.";
      throw std::logic_error(os.str());
    }

    this->R = par.size() / (k*k + k*k);

    this->Sigma = cube(&par[0], k, k, R);
    this->Sigmae = cube(&par[k*k*R], k, k, R);

    this->V = cube(k, k, this->ref_tree_.num_nodes());
    this->V_1 = cube(k, k, this->ref_tree_.num_nodes());
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

      V.slice(i) = ti * Sigma.slice(ri); 

      if(i < this->ref_tree_.num_tips()) {
        V.slice(i) += Sigmae.slice(ri);
      }

      arma::uvec ui(1);
      ui(0) = i;

      V_1.slice(i)(ki, ki) = inv(V.slice(i)(ki,ki));
      
      this->A.slice(i)(ki,ki) = -0.5*V_1.slice(i)(ki,ki);
      this->b(ki,ui).fill(0.0);
      C.slice(i)(kj,kj) = -0.5 * I(ki,kj).t() * V_1.slice(i)(ki,ki) * I(ki,kj);
      d(kj,ui).fill(0.0);
      E.slice(i)(kj,ki) = I(ki,kj).t() * V_1.slice(i)(ki,ki);
      f(i) = -0.5*(ki.n_elem * M_LN_2PI + log(det(V.slice(i)(ki,ki))));
    }
  }
};


typedef splittree::TraversalTask<BM> QuadraticPolynomialBM;
}

#endif // QuadraticPolynomial_BM_H_

/*
 *  QuadraticPolynomial.h
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
#ifndef PCMBase_QuadraticPolynomial_H_
#define PCMBase_QuadraticPolynomial_H_

#include "splittree.h"
#include <armadillo>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <sstream>

namespace PCMBaseCpp {

template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of big data objects

  // tip-names correponding to the columns in Pc_
  std::vector<NameType> const& names_;

  arma::mat const& X_;

  // column vectors of present coordinates (0: missing, !=0: present) for each
  // tip in the tree. The number of rows in the matrix corresponds to the number
  // of traits.
  arma::imat const& Pc_;
  bool internal_pc_full_;

  
  NumericTraitData(
    std::vector<NameType> const& names,
    arma::mat const& X,
    arma::imat const& Pc, 
    bool internal_pc_full): names_(names), X_(X), Pc_(Pc), internal_pc_full_(internal_pc_full) {}
};

struct LengthAndRegime {
  double length_;
  arma::uword regime_;

  LengthAndRegime() {}

  LengthAndRegime(double length, arma::uword regime):
    length_(length), regime_(regime) {}
};

struct LengthRegimeAndJump {
  double length_;
  arma::uword regime_;
  arma::u8 jump_;
  
  LengthRegimeAndJump() {}
  
  LengthRegimeAndJump(double length, arma::uword regime, arma::u8 jump):
    length_(length), regime_(regime), jump_(jump) {}
};

template<class RegimeType>
std::vector<arma::uword> mapRegimesToIndices(
    std::vector<RegimeType> const& regimes,
    std::vector<RegimeType> const& regimes_unique) {
  
  if(regimes_unique.size() == 0) {
    throw std::logic_error("ERR:03101:PCMBaseCpp:QuadraticPolynomial.h:mapRegimesToIndices:: regimes_unique has 0 length but should have at least one regime.");
  }
  std::unordered_map<RegimeType, arma::uword> map_regimes;
  arma::uword next_regime = 0;
  for(auto r: regimes_unique) {
    auto it = map_regimes.insert(std::pair<RegimeType, arma::uword>(r, next_regime));
    if(!it.second) {
      std::ostringstream os;
      os<<"ERR:03102:PCMBaseCpp:QuadraticPolynomial.h:mapRegimesToIndices:: The regime named '"<<r<<"' is dupliclated. Remove duplicates from regimes_unique.";
      throw std::logic_error(os.str());
    } else {
      ++next_regime;
    }
  }
  std::vector<arma::uword> regimeIndices;
  for(auto r: regimes) {
    auto it = map_regimes.find(r);
    if(it == map_regimes.end()) {
      std::ostringstream os;
      os<<"ERR:03103:PCMBaseCpp:QuadraticPolynomial.h:mapRegimesToIndices:: The regime named '"<<r<<"' was not found in regimes_unique.";
      throw std::logic_error(os.str());
    } else {
      regimeIndices.push_back(it->second);
    }
  }
  return regimeIndices;
}

template<class Tree>
class PresentCoordinates: public splittree::TraversalSpecification<Tree> {
public:
  typedef PresentCoordinates<Tree> MyType;
  typedef splittree::TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef splittree::PostOrderTraversal<MyType> AlgorithmType;
  typedef int ParameterType; // dummy ParameterType
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef arma::uvec StateType;

  arma::imat Pc_;
  bool internal_pc_full_;

  PresentCoordinates(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {
    if(input_data.Pc_.n_cols != this->ref_tree_.num_tips()) {
      throw std::invalid_argument(
          "ERR:03111:PCMBaseCpp:QuadraticPolynomial.h:PresentCoordinates:: The input matrix Pc_ must have as many rows as the number of traits and as many columns as the number of tips.");
    } else {
      // We need to be careful with the typedefs splittree::uvec and arma::uvec.
      using namespace arma;

      // number of tips
      uword M = this->ref_tree_.num_nodes();
      uword N = this->ref_tree_.num_tips();

      // number of traits
      uword k = input_data.Pc_.n_rows;

      if(k == 0) {
        throw std::invalid_argument(
            "ERR:03112:PCMBaseCpp:QuadraticPolynomial.h:PresentCoordinates:: The input matrix Pc_ must have as many rows as the number of traits. The number of traits should be at least 1 but was 0.");
      }

      this->internal_pc_full_ = input_data.internal_pc_full_;
      
      if(internal_pc_full_) {
        this->Pc_ = imat(k, M, fill::ones);
      } else {
        this->Pc_ = imat(k, M, fill::zeros);
      }

      uvec ordNodes(
          this->ref_tree_.OrderNodesPosType(
              input_data.names_, static_cast<uword>(splittree::NA_UINT)));

      Pc_.cols(0, N - 1) = input_data.Pc_.cols(ordNodes);
    }
  }

  void SetParameter(ParameterType const& par) {}

  void PruneNode(uint i, uint i_parent) {
    if(!internal_pc_full_) {
      using namespace arma;
      Pc_.col(i_parent) += Pc_.col(i);  
    }
  }

  // Present coordinates for a node (to be called after tree traversal)
  arma::uvec PcForNodeId(uint i) {
    using namespace arma;
    std::vector<arma::uword> pc;
    for(uword u = 0; u < Pc_.n_rows; ++u) {
      if(Pc_(u, i) > 0) {
        pc.push_back(u);
      }
    }
    return arma::uvec(&pc[0], pc.size());
  }

  StateType StateAtRoot() {
    return PcForNodeId(this->ref_tree_.num_nodes() - 1);
  }
};

template<class Tree>
class QuadraticPolynomial: public splittree::TraversalSpecification<Tree> {
public:
  typedef splittree::TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef std::vector<double> StateType;

  typedef splittree::TraversalTaskLightweight<
    PresentCoordinates<TreeType> > PresentCoordinatesTask;
  
  // singularity threshold for the determinant of V_i
  double threshold_SV_ = 1e-6;
  
  //
  // Input data consists of multiple trait values for each tip. Each
  // column corresponds to a tip.
  //
  arma::mat X;

  //
  // Coefficients used to calculate L, m, r, which are calculated in the
  // InitNode method of derived classes.
  //
  arma::cube A;
  arma::mat b;
  arma::cube C;
  arma::mat d;
  arma::cube E;
  arma::vec f;
  
  //
  // linear intercept and coefficient for the condiation mean: 
  // E[X_i|X_j] = omega + Phi X_j
  arma::mat omega;
  arma::cube Phi;
  
  // 
  // Variance - covariance matrices
  //
  arma::cube V;
  arma::cube V_1;
  

  //
  // Coefficients of the quadratic polynomial for each node
  // X.col(i).t() * L.slice(i) * X.col(i) + m.col(i) * X.col(i) + r(i)
  //
  arma::cube L;
  arma::mat m;
  arma::vec r;

  // present coordinates
  std::vector<arma::uvec> pc;

  // number of traits
  arma::uword k;

  QuadraticPolynomial(
    TreeType const& tree,
    NumericTraitData<typename TreeType::NodeType> const& input_data):

    BaseType(tree),

    X(input_data.X_),

    // all these fields have to be initialized with 0 during SetParameter.
    A(X.n_rows, X.n_rows, tree.num_nodes()),
    b(X.n_rows, tree.num_nodes()),
    C(X.n_rows, X.n_rows, tree.num_nodes()),
    d(X.n_rows, tree.num_nodes()),
    E(X.n_rows, X.n_rows, tree.num_nodes()),
    f(tree.num_nodes()),

    omega(X.n_rows, tree.num_nodes()),
    Phi(X.n_rows, X.n_rows, tree.num_nodes()),
    
    V(X.n_rows, X.n_rows, tree.num_nodes()),
    V_1(X.n_rows, X.n_rows, tree.num_nodes()),
  
    L(X.n_rows, X.n_rows, tree.num_nodes()),
    m(X.n_rows, tree.num_nodes()),
    r(tree.num_nodes()),

    k(X.n_rows) {

    arma::uvec ordNodes(
        this->ref_tree_.OrderNodesPosType(
            input_data.names_, static_cast<arma::uword>(splittree::NA_UINT)));

    this->X.cols(0, this->ref_tree_.num_tips() - 1) = X.cols(ordNodes);

    PresentCoordinatesTask pc_task(tree, input_data);
    pc_task.TraverseTree(0, 1);

    for(size_t i = 0; i < tree.num_nodes(); ++i) {
      pc.push_back(pc_task.spec().PcForNodeId(i));
      if(pc[i].n_elem == 0) {
        std::ostringstream oss;
        oss<<"ERR:03121:PCMBaseCpp:QuadraticPolynomial.h:QuadraticPolynomial:: Some tips ("<< this->ref_tree_.FindNodeWithId(i) <<") have 0 present coordinates. Consider removing these tips.";
        throw std::logic_error(oss.str());
      }
    }
  }

  StateType StateAtNode(arma::uword i) const {
    using namespace std;

    StateType res(k*k + k + 1 + k*k + k + k*k + k + k*k + 1 + k + k*k + k*k + k*k);

    copy(L.begin_slice(i), L.end_slice(i), res.begin());
    copy(m.begin_col(i), m.end_col(i), res.begin() + k*k);
    res[k*k + k] = r(i);
    copy(A.begin_slice(i), A.end_slice(i), res.begin() + k*k + k + 1);
    copy(b.begin_col(i), b.end_col(i), res.begin() + k*k + k + 1 + k*k);
    copy(C.begin_slice(i), C.end_slice(i), res.begin() + k*k + k + 1 + k*k + k);
    copy(d.begin_col(i), d.end_col(i), res.begin() + k*k + k + 1 + k*k + k + k*k);
    copy(E.begin_slice(i), E.end_slice(i), res.begin() + k*k + k + 1 + k*k + k + k*k + k);
    res[k*k + k + 1 + k*k + k + k*k + k + k*k] = f(i);
    copy(omega.begin_col(i), omega.end_col(i), res.begin() + k*k + k + 1 + k*k + k + k*k + k + k*k + 1);
    copy(Phi.begin_slice(i), Phi.end_slice(i), res.begin() + k*k + k + 1 + k*k + k + k*k + k + k*k + 1 + k);
    copy(V.begin_slice(i), V.end_slice(i), res.begin() + k*k + k + 1 + k*k + k + k*k + k + k*k + 1 + k + k*k);
    copy(V_1.begin_slice(i), V_1.end_slice(i), res.begin() + k*k + k + 1 + k*k + k + k*k + k + k*k + 1 + k + k*k + k*k);
    
    return res;
  }

  StateType StateAtRoot() const {
    return StateAtNode(this->ref_tree_.num_nodes() - 1);
  }

  // this function is to be called by daughter classes only after they have
  // initialized omega, Phi and V
  inline void CalculateAbCdEf(uint i) {
    using namespace arma;
    splittree::uint j = this->ref_tree_.FindIdOfParent(i);
    uvec kj = pc[j], ki = pc[i];
    uvec ui(1);
    ui(0) = i;
    
    A.slice(i)(ki,ki) = -0.5 * V_1.slice(i)(ki,ki);
    //A[ki,ki,i] <- (-0.5*V_1[ki,ki,i])
    E.slice(i)(kj,ki) = Phi.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki);
    //E[kj,ki,i] <- t(Phi[ki,kj,i]) %op% V_1[ki,ki,i]
    b(ki,ui) = V_1.slice(i)(ki,ki) * omega(ki,ui);
    //b[ki,i] <- V_1[ki,ki,i] %*% omega[ki,i]
    C.slice(i)(kj,kj) = -0.5 * E.slice(i)(kj,ki) * Phi.slice(i)(ki,kj);
    //C[kj,kj,i] <- -0.5 * matrix(E[kj,ki,i], sum(kj), sum(ki)) %*% matrix(Phi[ki,kj,i], sum(ki), sum(kj))
    d(kj,ui) = -E.slice(i)(kj,ki) * omega(ki,ui);
    //d[kj,i] <- -E[kj,ki,i] %op% omega[ki,i]
    f(i) = -0.5*(ki.n_elem * M_LN_2PI + log(det(V.slice(i)(ki,ki))) +
      omega(ki,ui).t() * V_1.slice(i)(ki,ki) * omega(ki,ui) ).at(0,0);
    //  (Theta.col(ri).t() * (I.rows(ki)-e_Ht.slice(i).rows(ki)).t() + xi*mj.col(ri).t()*e_Ht.slice(i).rows(ki).t()) * V_1.slice(i)(ki,ki) * ((I.rows(ki)-e_Ht.slice(i).rows(ki)) * Theta.col(ri) + xi*e_Ht.slice(i).rows(ki)*mj.col(ri))).at(0,0);
    //f[i] <- -0.5 * (t(omega[ki,i]) %*% V_1[ki,ki,i] %*% omega[ki,i] + sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i])))) 
  }
  
  inline void InitNode(uint i) {
    L.slice(i).fill(0.0);
    m.col(i).fill(0.0);
    r(i) = 0.0;
  }

  inline void VisitNode(uint i) {
    using namespace arma;
    using namespace std;

    splittree::uint j = this->ref_tree_.FindIdOfParent(i);

    uvec kj = pc[j], ki = pc[i];

    // check that V.slice(i)(ki,ki) is non-singular
    vec svd_V = svd(V.slice(i)(ki,ki));
    double ratio_SV = (*(svd_V.cend()-1))/(*svd_V.cbegin());
    if(ratio_SV < threshold_SV_) {
      ostringstream oss;
      oss<<"ERR:03131:PCMBaseCpp:QuadraticPolynomial.h:VisitNode:: The matrix V for node "<<
        this->ref_tree_.FindNodeWithId(i)<<" is nearly singular: min(svd_V)/max(svd_V)="<<
          ratio_SV<<", det(V)="<<det(V.slice(i)(ki,ki))<<
            ". Check the model parameters and the length of the branch"<<
              " leading to the node. For details on this error, read the User Guide.";
      throw logic_error(oss.str());
    }
    
    if(i < this->ref_tree_.num_tips()) {
      arma::uvec ui(1);
      ui(0) = i;

      // ensure symmetry of L.slice(i)
      L.slice(i) = 0.5 * (C.slice(i) + C.slice(i).t());
      r(i) = (X(ki,ui).t() * A.slice(i)(ki,ki) * X(ki,ui) +
        X(ki,ui).t() * b(ki,ui) + f(i)).at(0,0);
      m(kj, ui) = d(kj, ui) + E.slice(i)(kj,ki) * X(ki,ui);
    } else {
      uvec ui(1);
      ui(0) = i;

      mat AplusL = A.slice(i)(ki,ki) + L.slice(i)(ki,ki);
      // Ensure symmetry of AplusL:
      AplusL = 0.5 * (AplusL + AplusL.t());

      mat AplusL_1 = inv(AplusL);
      mat EAplusL_1 = E.slice(i)(kj,ki) * AplusL_1;
      double logDetVNode = log(det(-2*AplusL));

      r(i) = (f(i) + r(i) + 0.5 * ki.n_elem * M_LN_2PI - 0.5 * logDetVNode -
        0.25 * (b(ki,ui) + m(ki,ui)).t() * AplusL_1 * (b(ki,ui) + m(ki,ui))).at(0,0);
      m(kj,ui) = d(kj,ui) - 0.5*EAplusL_1 * (b(ki,ui) + m(ki,ui));
      L.slice(i)(kj,kj) = C.slice(i)(kj,kj) - 0.25 * EAplusL_1 * E.slice(i)(kj,ki).t();
      // Ensure symmetry of L.slice(i)
      L.slice(i)(kj,kj) = 0.5 * (L.slice(i)(kj,kj) + L.slice(i)(kj,kj).t());
    }
  }

  inline void PruneNode(uint i, uint i_parent) {
    L.slice(i_parent) += L.slice(i);
    m.col(i_parent) += m.col(i);
    r(i_parent) += r(i);
  }
};
}
#endif // PCMBase_QuadraticPolynomial_H_

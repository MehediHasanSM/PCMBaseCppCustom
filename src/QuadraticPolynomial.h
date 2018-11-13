/*
 *  QuadraticPolynomial.h
 *  PCMBaseCpp
 *
 * Copyright 2017,2018 Venelin Mitov
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

#include "SPLITT.h"
#include <armadillo>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <iostream>


namespace PCMBaseCpp {

// Check if a square matrix is diagonal
template<class MatType>
inline bool IsDiagonal(MatType const& X) {
  using namespace arma;
  
  for(int i = 0; i < X.n_rows; ++i) {
    for(int j = i + 1; j < X.n_cols; ++j) {
      if( X(i, j) != 0.0 || X(j, i) != 0.0 ) {
        return false;
      }
    }
  }
  return true;
}

template<class MatType>
inline bool IsSingular(MatType const& X, double threshold_SV) {
  using namespace arma;
  vec svd_V = svd(X);
  double ratio_SV = (*(svd_V.cend()-1)) / (*svd_V.cbegin());
  return (!std::isfinite(ratio_SV) || ratio_SV < threshold_SV);
}

template<class MatType, class VecType>
inline void PairSums(MatType& pairSums, VecType const& elems) {
  using namespace arma;
  uword k = elems.n_elem;
  for(uword i = 0; i < k; ++i)
    for(uword j = i; j < k; ++j) 
      pairSums(i,j) = pairSums(j,i) = elems(i) + elems(j);
}

template<class MatEigvalType, class CubeEigvecType, class CubeHType>
inline void DecomposeH(MatEigvalType& lambda, CubeEigvecType& P, CubeEigvecType& P_1, CubeHType const& H, 
                       arma::uword r, double threshold_SV) {
  using namespace arma;
  cx_vec eigval(H.slice(r).n_rows);
  eigval.fill(std::complex<double>(0.0, 0.0));
  cx_mat eigvec(H.slice(r).n_rows, H.slice(r).n_rows);
  
  if(IsDiagonal(H.slice(r))) {
    // This is a workaround that was needed specifically for scalar-diagonal matrices H and Sigma_x
    eigval.set_real(real(H.slice(r).diag()));
    eigvec = eye<cx_mat>(H.slice(r).n_rows, H.slice(r).n_rows);
  } else {
    eig_gen(eigval, eigvec, H.slice(r));
  }

  lambda.col(r) = eigval;
  P.slice(r) = eigvec;
  if(IsSingular(P.slice(r), threshold_SV)) {
    std::ostringstream os;
    os<<"ERR:03402:PCMBaseCpp:QuadraticPolynomial.h:DecomposeH:: Defective H matrix:"<<
      H.slice(r)<<" - the matrix of eigenvectors is computationally singular.";
    throw std::logic_error(os.str());
  }
  P_1.slice(r) = inv(P.slice(r));
  
}


template<class MatType>
inline void CDFExpDivLambda(MatType& fLambda_ij, MatType const& Lambda_ij, double time, double threshold_Lambda_ij) {
  using namespace arma;
  for(uword w = 0; w < Lambda_ij.n_cols; ++w) {
    for(uword v = w; v < Lambda_ij.n_cols; ++v) {
      if(abs(Lambda_ij(w,v)) < threshold_Lambda_ij) {
        fLambda_ij(w,v) = fLambda_ij(v,w) = time;
      } else {
        fLambda_ij(w,v) = fLambda_ij(v,w) =
          (1.0 - exp(-Lambda_ij(w,v) * time)) / Lambda_ij(w,v);
      }
    }
  }  
}

template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of big data objects

  // tip-names correponding to the columns in Pc_
  std::vector<NameType> const& names_;

  arma::mat const& X_;

  // column vectors of present coordinates (0: missing, !=0: present) for each
  // node in the tree. 
  std::vector<arma::uvec> const& Pc_;
  
  uint k_;
  uint R_;
  
  std::vector<std::string> regime_models_;
  double threshold_SV_;
  double threshold_EV_;
  double threshold_skip_singular_;
  double threshold_Lambda_ij_;
  bool skip_singular_;
  
  NumericTraitData(
    std::vector<NameType> const& names,
    arma::mat const& X,
    std::vector<arma::uvec> const& Pc, 
    uint R,
    std::vector<std::string> regime_models,
    double threshold_SV,
    double threshold_EV,
    double threshold_skip_singular,
    bool skip_singular,
    double threshold_Lambda_ij): names_(names), X_(X), Pc_(Pc), k_(X.n_rows), 
      R_(R), regime_models_(regime_models),
      threshold_SV_(threshold_SV), 
      threshold_EV_(threshold_EV), 
      threshold_skip_singular_(threshold_skip_singular),
      threshold_Lambda_ij_(threshold_Lambda_ij),
      skip_singular_(skip_singular) {}
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

// Conditional Gaussian distribution of trait vector at a daughter 
// node, Xi, given trait vector at its parent, Xj, assuming that the conditional mean
// depends linearly in Xj, i.e. Mean(Xi) = omega + Phi %*% Xj, and the variance does
// not depend on Xj. This is an abstract class implemented by specific models
class CondGaussianOmegaPhiV {
public:
  virtual arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) = 0;
  virtual void CalculateOmegaPhiV(uint i, arma::uword ri, arma::mat& omega, arma::cube& Phi, arma::cube& V) = 0;
  virtual ~CondGaussianOmegaPhiV() {}
};

template<class Tree>
class QuadraticPolynomial: public SPLITT::TraversalSpecification<Tree> {
public:
  typedef SPLITT::TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef std::vector<double> StateType;

  // singular value threshold for the determinant of V_i
  double threshold_SV_ = 1e-6;
  // positive eigenvalue threshold for V_i. 
  double threshold_EV_ = 1e-4;
  
  // threshold specifying the maximum allowed branch length for skipping the branch
  // if the corresponding matrix V is singular. This option matters only if
  // skip_singular_ is true.
  double threshold_skip_singular_ = 1e-4;
  
  // denotes branches for which V is singular
  // the treatment of these branches depends on the option PCMBase.Singular.Skip
  // If this option is set to TRUE (default), then these branches are treated as
  // 0-length and the L,m,r values accumulated from their children are added 
  // up to their parent branches without modification.
  // ATTENTION: using std::vector<bool> instead of std::vector<int> is causing a 
  // bug in parallel mode (omp for simd). This was found in the unit test for 
  // the White model, test-White.R in PCMBase, using the intel compiler icpc 
  // (icpc version 17.0.5 (gcc version 4.9.0 compatibility)), command line:
  // icpc -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include -I/usr/local/include/freetype2 -I/opt/X11/include -I"/Users/vmitov/Library/R/3.3/library/Rcpp/include" -I"/Users/vmitov/Library/R/3.3/library/RcppArmadillo/include"   -fPIC  -Wall -mtune=core2 -g -O2  -std=c++11 -fopenmp -Wall -O2 -march=native -c Rcpp.cpp -o Rcpp.o
  // the problem appears to be solved when using std::vector<int>.
  std::vector<int> singular_branch_;
  
  bool skip_singular_;
  
  
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

  std::vector<CondGaussianOmegaPhiV*> ptr_cond_dist_;
  
  QuadraticPolynomial(
    TreeType const& tree,
    NumericTraitData<typename TreeType::NodeType> const& input_data):

    BaseType(tree),

    threshold_SV_(input_data.threshold_SV_),
    threshold_EV_(input_data.threshold_EV_),
    threshold_skip_singular_(input_data.threshold_skip_singular_),
    singular_branch_(tree.num_nodes(), false),
    skip_singular_(input_data.skip_singular_),
    
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

    arma::uvec ordTips(
        this->ref_tree_.OrderNodesPosType(
            input_data.names_, static_cast<arma::uword>(SPLITT::G_NA_UINT)));

    this->X.cols(0, this->ref_tree_.num_tips() - 1) = X.cols(ordTips);

    //PresentCoordinatesTask pc_task(tree, input_data);
    //pc_task.TraverseTree(0, 1);

    SPLITT::uvec node_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), tree.num_nodes());
    arma::uvec ordNodes(
        this->ref_tree_.OrderNodesPosType(
            node_names, static_cast<arma::uword>(SPLITT::G_NA_UINT)));
    
    for(size_t i = 0; i < tree.num_nodes(); ++i) {
      //pc.push_back(pc_task.spec().PcForNodeId(i));
      pc.push_back(input_data.Pc_[ordNodes(i)]);
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
    SPLITT::uint j = this->ref_tree_.FindIdOfParent(i);
    uvec kj = pc[j], ki = pc[i];
    uvec ui(1);
    ui(0) = i;
    
    A.slice(i)(ki,ki) = -0.5 * V_1.slice(i)(ki,ki);
    E.slice(i)(kj,ki) = Phi.slice(i)(ki,kj).t() * V_1.slice(i)(ki,ki);
    b(ki,ui) = V_1.slice(i)(ki,ki) * omega(ki,ui);
    C.slice(i)(kj,kj) = -0.5 * E.slice(i)(kj,ki) * Phi.slice(i)(ki,kj);
    d(kj,ui) = -E.slice(i)(kj,ki) * omega(ki,ui);
    f(i) = -0.5*(ki.n_elem * M_LN_2PI + log(det(V.slice(i)(ki,ki))) +
      omega(ki,ui).t() * V_1.slice(i)(ki,ki) * omega(ki,ui) ).at(0,0);
  }
  
  inline void InitLmr(uint i) {
    L.slice(i).fill(0.0);
    m.col(i).fill(0.0);
    r(i) = 0.0;
    singular_branch_[i] = 0;
  }
  
  inline void InitNode(uint i) {
    using namespace arma;
    using namespace std;
    InitLmr(i);
    
    if(i < this->ref_tree_.num_nodes() - 1) {
      
      auto ri = this->ref_tree_.LengthOfBranch(i).regime_;
      auto ti = this->ref_tree_.LengthOfBranch(i).length_;
      
      if(ptr_cond_dist_.size() == 1) {
        ptr_cond_dist_[0]->CalculateOmegaPhiV(i, ri, omega, Phi, V);
      } else {
        ptr_cond_dist_[ri]->CalculateOmegaPhiV(i, 0, omega, Phi, V);
      }
      
      arma::uvec ki = pc[i];
      
      // force symmetry of V, which could be lost due to numerical inprecision
      V.slice(i)(ki,ki) = 0.5 * (V.slice(i)(ki,ki) + V.slice(i)(ki,ki).t());
      
      if( IsSingular(V.slice(i)(ki,ki), threshold_SV_) ) {
        singular_branch_[i] = 1;
        if(!skip_singular_ || ti > threshold_skip_singular_) {
          ostringstream oss;
          oss<<"ERR:03131:PCMBaseCpp:QuadraticPolynomial.h:InitNode:: The matrix V for node "<<
            this->ref_tree_.FindNodeWithId(i)<<" is nearly singular: "<<V.slice(i)(ki,ki)<<
                ". Check the model parameters and the length of the branch"<<
                  " leading to the node. For details on this error, read the User Guide.";
          throw logic_error(oss.str());  
        } 
      } 
      
      if(!singular_branch_[i]) {
        // Check V is positive definite: all eigen-values must be strictly positive
        cx_vec eigval(ki.n_elem);
        eigval.fill(std::complex<double>(0.0, 0.0));
        cx_mat eigvec(ki.n_elem, ki.n_elem);
        
        arma::mat VSliceIki = V.slice(i)(ki,ki);
        if(IsDiagonal(VSliceIki)) {
          // this is a workaround needed apparently only for diagonal matrices
          // with duplicated elements on the diagonal, but should work fine 
          // with any diagonal matrix
          eigval.set_real(VSliceIki.diag());
          eigvec = eye<cx_mat>(ki.n_elem, ki.n_elem);
        } else {
          eig_gen(eigval, eigvec, VSliceIki);
        }

        
        // eig_gen(eigval, eigvec, V.slice(i)(ki,ki));
        vec re_eigval = real(eigval);
        
        for(double eigv: re_eigval) {
          if(eigv < threshold_EV_) {
            ostringstream oss;
            oss<<"ERR:03132:PCMBaseCpp:QuadraticPolynomial.h:InitNode:: The matrix V for node "<<
              this->ref_tree_.FindNodeWithId(i)<<
                " is nearly singular or not positive definite; near-0 or negative eigenvalue found: "<<eigv<<
                "V.slice(i)(ki,ki): "<<V.slice(i)(ki,ki)<<". Check the model parameters.";
            throw logic_error(oss.str());
          }
        }
        
        V_1.slice(i)(ki, ki) = real(eigvec * diagmat(1/eigval) * eigvec.t());
        //V_1.slice(i)(ki, ki) = inv(V.slice(i)(ki,ki));
        //V_1.slice(i)(ki, ki) = inv_sympd(V.slice(i)(ki,ki));
        CalculateAbCdEf(i);  
      }
    }  
    
  }

  inline void VisitNode(uint i) {
    using namespace arma;
    using namespace std;

    if(!singular_branch_[i]) {
      SPLITT::uint j = this->ref_tree_.FindIdOfParent(i);
  
      uvec kj = pc[j], ki = pc[i];
  
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
  }

  inline void PruneNode(uint i, uint i_parent) {
    L.slice(i_parent) += L.slice(i);
    m.col(i_parent) += m.col(i);
    r(i_parent) += r(i);
  }
};
}

#endif // PCMBase_QuadraticPolynomial_H_
